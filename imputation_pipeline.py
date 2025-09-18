import argparse
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import gzip
import glob
import numpy as np
import os
import pandas as pd
from pathlib import Path
import re
import shutil
import subprocess

# last updated Sep. 14, 2025
# extended from initial in-house implementation of TOPMed pipeline by Minhui Chen
# This script will also use bcftools and tabix!
# bcftools and tabix can be loaded in slurm script via "module load bcftools" and "module load htslib"

# TODO
    #Add in debug flag logic


#
# Program paths
#
PLINK = '/project/haiman_625/Software/imputation_pipeline_CharlestonGroup/programs/plink'
EAGLE = '/project/haiman_625/Software/imputation_pipeline_CharlestonGroup/programs/eagle'
GMAP = '/project/haiman_625/Software/imputation_pipeline_CharlestonGroup/programs/genetic_map_hg38_withX.txt.gz'
MM4 = '/project/chia657_28/programs/minimac4/bin/minimac4'
HIST_R_SCRIPT = 'script_to_plot_imputed_variants_vs_ref.R'


parser = argparse.ArgumentParser()
parser.add_argument('--tar', required=True, help='target data vcf (must be compressed)')
parser.add_argument('--ref', required=True, help='phased reference VCF (must be compressed)')
parser.add_argument('--refm3', required=True, help='reference m3vcf')
parser.add_argument('--chrom', required=True, help='format the same way as hg38 (e.g. "chr22")')
parser.add_argument('--outdir', required=True, help='directory to write to (without trailing "/")')
parser.add_argument('--workers', required=True, help='number of cpus to use in parallel (e.g. num. chunks to work on at same time)')
parser.add_argument('--mem', required=True, help='mem in gb to use')
parser.add_argument('--windowsize', required=False, help='mm4 window size in bp', default=500000)
parser.add_argument('--chunksize', required=False, help='chunk size in bp', default=25e6)
parser.add_argument('--overlap', required=False, help='overlap size in bp', default=5e6)
parser.add_argument('--multicore', required=False, help='number of workers to use when multithreading (e.g. number of cpus to use during imputation)', default=5)
parser.add_argument('--nophasing', default=False, action='store_true', help='Skips phasing during pipeline (for analyses such as meta-imputation)')
parser.add_argument('--debug', default=False, action='store_true', help='Enable debugging: skips running phasing/imputation subprocesses but still writes log files.')
parser.add_argument('--noplot', default=False, action='store_true', help='Skip plotting of histogram of imputed vs. reference variants.')

#
# Parse args
#

args = parser.parse_args()
tar = args.tar
ref = args.ref
refm3 = args.refm3
chrom = args.chrom
outdir = args.outdir
workers = int(args.workers)
mem = float(args.mem)
chunksize = int(args.chunksize)
overlap = int(args.overlap)
windowsize = int(args.windowsize)
multicore = int(args.multicore)
nophasing = args.nophasing
debug_mode = args.debug
skiphistplot = args.noplot

#directories
tmpdir = outdir+f'/tmp_{chrom}' #to be deleted after successful imputation!
qcdir = tmpdir+'/qc'
chr_dir = f'{qcdir}/{chrom}'
impute_chr_dir = tmpdir+f'/imputation/{chrom}'

#files
posfp = tmpdir+f'/{chrom}.pos.txt'

#
# Params
#

mem_per_worker = int((mem/workers)*1000)
max_multithread = max(1,int(workers/multicore))
min_valid_variants = 3
min_percent_valid_variants = 0.5 #check vs reference, min 50%
minsampleCallRate = 0.5
alleles = ['A', 'T', 'C', 'G']
minsnpCallRate = 0.9
min_maf = 1e-10 # remove monomorphic snps
buffersize = int(overlap/2)
min_window_offset = 1000000  #shift 1Mb at the beginning of chr


# Write/overwite log for run
log = f'{outdir}/{chrom}.log.txt'


#
# Functions
#


def write_log(fp,string):
    with open(fp, 'w') as loghandle:
        loghandle.write(string)


def append_log(fp,string):
    with open(fp, 'a') as loghandle:
        loghandle.write(string)


def position_windows(pos, size, start=None, stop=None, step=None):
    """
    Break up variants into windows from pos
    Window is shifted if no variants are found in the beginning of the window
    Added 'final' chunk that has snps spanning entire interval
        updated trimming logic to calculate midpoint of each chunk on-the-fly
        chunk is expected to impute better where variants don't span entire window
    """

    #final chunk
    finalwindow_ceiling = pos[-1]//min_window_offset + 1
    final_stop = finalwindow_ceiling*min_window_offset
    final_start = final_stop - chunksize + 1

    # determine start and stop positions
    if start is None:
        start = pos[0]
    if stop is None:
        stop = final_stop
    if step is None:
        # non-overlapping
        step = size

    windows = []
    window_start = start
    window_stop = start + size - 1

    while (window_start < final_start):
        snps_in_beginning = False
        num_snps_in_first_section = [x for x in pos if window_start <= x <= window_start + min_window_offset - 1]
        if len(num_snps_in_first_section) > 0:
            snps_in_beginning = True

        #increment logic
        if snps_in_beginning:
            windows.append([window_start, window_stop])
            window_start += step
            window_stop += step
        else:
            #shift window
            append_log(log, f'{window_start}-{window_stop}: no snps detected from {window_start} to {window_start+min_window_offset} - shifted window {min_window_offset}\n')
            window_start += min_window_offset
            window_stop += min_window_offset

        #early stop condition to avoid 2nd-to-last chunk highly overlapping last chunk
        if final_start - window_start < overlap:
            break

    windows.append([final_start, final_stop])
    return np.asarray(windows)


def get_chrompos_pd_from_vcf(filepath):
    #read in vcf and parse for chrom, pos, ref, and alt

    lines = [l.split() for l in subprocess.check_output(f"zcat {filepath} | grep ^[^#] | cut -f1-2,4-5", shell=True).decode().splitlines()]
    all_pos = [l[1] for l in lines]
    lines = [f'{l[0]},{l[1]},{sorted([l[2],l[3]])[0]}{sorted([l[2],l[3]])[1]}' for l in lines]
    df = pd.DataFrame(lines, columns=['snp'])
    return df, all_pos


def create_chunks(window):
    #parse vcf between start and end
    start, stop = window
    chunkprefix = chr_dir+f'/chunk_{start}_{stop}'
    append_log(log,f'starting create_chunks {window}\n')
    result = subprocess.call(f'bcftools view -Oz -o {chunkprefix}.vcf.gz -r {chrom}:{start}-{stop} {tar}', shell=True)
    append_log(log,f'finished create_chunks {window}\n')
    return result


def create_pathdict_by_window(window):
    # Each chunk's output files are named corresponding to the window interval.

    start, stop = window

    d = {}

    d['chunkprefix'] = chr_dir+f'/chunk_{start}_{stop}'
    d['chunkvcf'] = chr_dir+f'/chunk_{start}_{stop}.vcf.gz'
    d['chrmissing_dir'] = f'{qcdir}/{chrom}_missing'

    d['chunkmissingprefix'] = f'{d["chrmissing_dir"]}/chunk_{start}_{stop}'
    d['chrfilter_dir'] = f'{qcdir}/{chrom}_filteredChunk'
    d['chrfilterchunk'] = f'{d["chrfilter_dir"]}/chunk_{start}_{stop}.txt'

    d['snpsfile'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.snps'

    d['bcf'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.bcf'
    d['index'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.bcf.csi'

    d['chrmissing_dir'] = f'{qcdir}/{chrom}_missing'
    d['chunkmissingprefix'] = f'{d["chrmissing_dir"]}/chunk_{start}_{stop}'

    d['phasedprefix'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.phased'
    d['phasedvcf'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.phased.vcf.gz'
    d['skipphasevcf'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.notrephased.vcf.gz'

    d['imputeddosefile'] = impute_chr_dir+f'/chunk_{start}_{stop}.imputed.dose.vcf.gz'
    d['imputedprefix'] = impute_chr_dir+f'/chunk_{start}_{stop}.imputed'

    return d


def qc_chunks(window):
    append_log(log,f'starting qc_chunks {window}\n')
    w_dict = create_pathdict_by_window(window)

    chunk_pd, _ = get_chrompos_pd_from_vcf(w_dict['chunkvcf'])
    num_valid_variants = ref_pd.loc[ref_pd['snp'].isin(chunk_pd['snp'])].shape[0]
    percent_valid_variants = num_valid_variants/ref_pd.shape[0]
    Path(w_dict['chrmissing_dir']).mkdir(parents=True, exist_ok=True)
    result = subprocess.call(f'{PLINK} --vcf {w_dict["chunkvcf"]} \
                                       --missing \
                                       --memory {mem_per_worker} \
                                       --out {w_dict["chunkmissingprefix"]}', shell=True, stdout=subprocess.DEVNULL)
    missing = pd.read_csv(w_dict['chunkmissingprefix']+'.imiss', sep='\s+')
    num_high_missing = missing.loc[missing['F_MISS'] > (1-minsampleCallRate)].shape[0]
    Path(w_dict['chrfilter_dir']).mkdir(parents=True, exist_ok=True)

    with open(w_dict['chrfilterchunk'], 'w') as f:
        if num_valid_variants < min_valid_variants or num_high_missing > 0:
            f.write(f'Exclude\t{num_valid_variants}\t{percent_valid_variants}\t{num_high_missing}')
        else:
            f.write('Keep')
    append_log(log,f'finished qc_chunks {window}\n')
    return result


def filter_chunks(window):
    w_dict = create_pathdict_by_window(window)

    #read bim into df
    lines = [l.split() for l in subprocess.check_output(f"zcat {w_dict['chunkvcf']} | grep ^[^#] | cut -f1-5", shell=True).decode().splitlines()]
    append_log(log,f'finished read chunkvcf into list {window}\n')

    lines = [[l[0], l[1], l[2], sorted([l[3],l[4]])[0], sorted([l[3],l[4]])[1]] for l in lines]
    bim = pd.DataFrame(lines, columns=['chr', 'pos', 'origsnpid','a1', 'a2'])
    append_log(log,f'finished create bim df from list {window}\n')
    append_log(log,f'window {window} datatypes: {bim.dtypes}\n')
    bim['pos'] = bim['pos'].astype(int)
    bim['snp'] = bim.apply(lambda row: f"{row['chr']},{row['pos']},{row['a1']}{row['a2']}", axis=1)
    append_log(log,f'finished create chrompos col in df {window}\n')

    # pick loc of snps only
    bim = bim.loc[(bim['a1'].isin(alleles+['0'])) & (bim['a2'].isin(alleles+['0']))]
    pos = np.array(bim['pos'])
    append_log(log,f'finished subset pos in df {window}\n')

    # make sure dis is not negative
    dis = pos[1:] - pos[:-1]
    if (dis < 0).sum() > 0:
        append_log(log,f'Positions are not inc in bim window: {window}\n')
        #return 0
        return 1

    snp = np.array(bim['snp'])
    duplicates = np.append(snp[:-1][dis<=1], snp[1:][dis<=1])
    bim = bim.loc[~bim['snp'].isin(duplicates)]

    append_log(log,f'finished remove dup pos in df {window}\n')

    # filter missing rate / call rate
    lmissfile = w_dict['chunkmissingprefix']+'.lmiss'
    highmissing = pd.read_csv(lmissfile, sep='\s+')
    append_log(log,f'finished parsing lmissfile {window}\n')
    highmissing = highmissing.loc[highmissing['F_MISS'] > (1-minsnpCallRate)]
    bim = bim.loc[~bim['snp'].isin(highmissing['SNP'])]
    append_log(log,f'finished remove highmissing in df {window}\n')

    # extract snps and convert to bcf, then index
    bim['origsnpid'].to_csv(w_dict['snpsfile'], sep='\t', index=False, header=False)
    append_log(log,w_dict['snpsfile']+'\n')
    bcftools_cmd = f'bcftools view --include ID==@{w_dict["snpsfile"]} -Ob -o {w_dict["bcf"]} {w_dict["chunkvcf"]}'
    result1 = subprocess.call(bcftools_cmd, shell=True)
    append_log(log,f'finished bcftools vcf to bcf {window}\n')
    result2 = subprocess.call(f'bcftools index {w_dict["bcf"]}', shell=True)
    append_log(log,f'finished bcftools indexing {window}\n')
    return result1 & result2


def phase_chunks(window):
    w_dict = create_pathdict_by_window(window)
    start, stop = window

    filteredchunk = open(w_dict['chrfilterchunk']).read().strip()
    if filteredchunk == 'Keep':
        cmd = f'{EAGLE} --vcfRef {ref} \
                                  --vcfTarget {w_dict["bcf"]} \
                                  --geneticMapFile {GMAP} \
                                  --chrom {chrom} \
                                  --bpStart {start} \
                                  --bpEnd {stop} \
                                  --outPrefix {w_dict["phasedprefix"]} \
                                  --numThreads {multicore} \
                                  --allowRefAltSwap \
                                  --vcfOutFormat z'

        cmd_for_log = "    \n".join(cmd.split())
        append_log(log,f'PHASING CHUNK {window}:\n{cmd_for_log}\n')
        result = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    else:
        result = subprocess.call(f'touch {w_dict["phasedvcf"]}', shell=True)
    return result


def convert_bcf_to_vcf_without_phasing(window):
    # previously plink -> bcf -> vcf
    # now vcf -> bcf -> vcf

    w_dict = create_pathdict_by_window(window)

    filteredchunk = open(w_dict['chrfilterchunk']).read().strip()
    if filteredchunk == 'Keep':
        result = subprocess.call(f'bcftools view -Oz -o {w_dict["skipphasevcf"]} --threads {multicore} {w_dict["bcf"]}', shell=True, stdout=subprocess.DEVNULL)
    else:
        result = subprocess.call(f'touch {w_dict["skipphasevcf"]}', shell=True)
    return result


def impute_chunks(window):
    w_dict = create_pathdict_by_window(window)
    start, stop = window

    if nophasing:
        preimputationvcf = w_dict['skipphasevcf']
    else:
        preimputationvcf = w_dict['phasedvcf']

    filteredchunk = open(w_dict['chrfilterchunk']).read().strip()
    append_log(log,f'impute chunks on {window}: {filteredchunk}\n')

    if filteredchunk == 'Keep' and stop - buffersize > start + buffersize:
        empty_imputation_region = not chunk_has_variants_to_parse(preimputationvcf, start, stop)

        if empty_imputation_region:
            cmd = f'touch {w_dict["imputeddosefile"]}'
            append_log(log,f'empty impute region: no variants between {start + buffersize} and {stop - buffersize}, running:{cmd}\n')
            result = subprocess.call(cmd, shell=True)
        else:
            cmd = f'{MM4} --refHaps {refm3} \
                                    --haps {preimputationvcf} \
                                    --chr {chrom} \
                                    --start {start} \
                                    --end {stop} \
                                    --window {windowsize} \
                                    --cpus {multicore} \
                                    --format GT,DS,GP \
                                    --ChunkLengthMb {int(chunksize/1000000)} \
                                    --ChunkOverlapMb {int(overlap/1000000)} \
                                    --allTypedSites \
                                    --meta \
                                    --noPhoneHome \
                                    --minRatio 0.00001 \
                                    --prefix {w_dict["imputedprefix"]}'
            cmd_for_log = "    \n".join(cmd.split())
            append_log(log,f'IMPUTING CHUNK {window}:\n{cmd_for_log}\n')
            result1 = subprocess.call(cmd, shell=True)
            if os.path.exists(f'{w_dict["imputeddosefile"]}.tbi'):
                os.remove(f'{w_dict["imputeddosefile"]}.tbi')
            result2 = subprocess.call(f'tabix -f -p vcf {w_dict["imputeddosefile"]}', shell=True, stdout=subprocess.DEVNULL)
            if os.path.exists(f'{w_dict["imputeddosefile"].replace("dose","empiricalDose")}.tbi'):
                os.remove(f'{w_dict["imputeddosefile"].replace("dose","empiricalDose")}.tbi')
            result3 = subprocess.call(f'tabix -f -p vcf {w_dict["imputeddosefile"].replace("dose","empiricalDose")}', shell=True, stdout=subprocess.DEVNULL)
            result = result1 and result2 and result3
    else:
        cmd = f'touch {w_dict["imputeddosefile"]}'
        append_log(log,f'chunk {window} failed "Keep" and start/stop check in impute chunk step, running:{cmd}\n')
        result = subprocess.call(cmd, shell=True)
    return result


def parse_info_between_pos(info_path, start_bp, end_bp):
    #read in info file
    #returns lines between the interval (inclusive, skips header) based on pos in snp name

    lines = []
    with open(info_path) as f:
        for rline in f.read().splitlines():
            if rline.startswith('SNP'):
                continue

            line = rline.split()
            pos = int(line[0].split(':')[1])
            if start_bp <= pos <= end_bp:
                lines.append(rline.strip())
    return lines


def chunk_has_variants_to_parse(fp, window_start, window_stop):
    #read lines of vcf to make sure there are enough variants that pass

    has_variants_after_trim = False

    for line in gzip.open(fp, 'rt'):
        if not re.search('^#', line):
            line = line.split('\t', 2)
            pos = int(line[1])

            if pos > window_start + buffersize and pos < window_stop - buffersize:
                has_variants_after_trim = True
                break

    return has_variants_after_trim


def process_chunk_windows_to_final_trim_windows(paths):
    #get start/stop of each window and save into list of tuples
    #hand off list to corresponding functions: vcf, empiricalDose, info

    path_windows = []
    for chunk in paths:
        cur_start = int(os.path.basename(chunk).split('_')[1])
        cur_stop = int(os.path.basename(chunk).split('_')[2].split('.')[0])
        path_windows.append([cur_start, cur_stop])

    max_idx = len(path_windows)-1
    idx = 0
    while idx < max_idx:
        cur_start, cur_stop = path_windows[idx]
        next_start, next_stop = path_windows[idx+1]

        if next_start < cur_stop:
            #overlap! Shift edges of both windows
            cur_overlap_midpoint = int((cur_stop+1-next_start)/2)
            new_cur_stop = cur_stop - cur_overlap_midpoint
            path_windows[idx] = [cur_start, new_cur_stop]
            new_next_start = next_start + cur_overlap_midpoint
            path_windows[idx+1] = [new_next_start, next_stop]

        idx += 1

    assert len(paths) == len(path_windows), 'Number of windows (%s) did not match number of paths given (%s)' % (len(path_windows), len(paths))
    return path_windows


def write_final_info_file(list_of_info, list_of_windows, output_file_name):
    #Takes in list of info filepaths and their corresponding ranges
    #Writes the subset of info lines to the given output filepath

    with gzip.open(output_file_name, 'wt') as g:
        g.write('SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1\n')
        for info_fp,window in zip(list_of_info,list_of_windows):
            cur_start, cur_stop = window
            valid_info_lines = parse_info_between_pos(info_fp, cur_start, cur_stop)
            for line in valid_info_lines:
                g.write(f'{line}\n')


def HELPER_write_final_vcf_file_tabix_and_trim(t):
    #parallelizable function to process vcf for concatenation
    #indexes original vcf, extract subset, and index subset
    #t is a tuple holding the original fp and region to be kept

    fp, start, stop = t
    new_fp = fp.replace('.vcf.gz','.finalsubset.vcf.gz')
    result1 = subprocess.call(f'tabix -f -p vcf {fp}', shell=True)
    result2 = subprocess.call(f'bcftools view -r {chrom}:{start}-{stop} -Oz -o {new_fp} {fp}', shell=True)
    result3 = subprocess.call(f'tabix -f -p vcf {new_fp}', shell=True)

    return result1 & result2 & result3


def write_final_vcf_file(list_of_vcfs, list_of_windows, output_file_name):
    #take in original list of vcfs, the ranges to trim them to, and write to output file
    #create list with new paths

    vcf_window_tuples = [[fp]+window for fp,window in zip(list_of_vcfs, list_of_windows)]

    append_log(log,f'Begin trimming: {list_of_vcfs}\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        trim_vcf_results = executor.map(HELPER_write_final_vcf_file_tabix_and_trim, vcf_window_tuples)
        executor.shutdown()

    for fp, result, window in zip(list_of_vcfs, trim_vcf_results, list_of_windows):
        append_log(log,f'Trimming {fp} to {window}: {result}\n')

    new_list_of_vcfs = [fp.replace('.vcf.gz','.finalsubset.vcf.gz') for fp in list_of_vcfs]
    chunkstring = ' '.join(new_list_of_vcfs)
    append_log(log,f'Begin concatenation of chunks: {chunkstring}\n')
    result1 = subprocess.call(f'bcftools concat --threads {workers} -Oz -o {output_file_name} -a -d all {chunkstring}', shell=True)
    append_log(log, f'Result of writing file {output_file_name}: {result1}\n')

    return result1


def HELPER_tabix_vcf_parallel(fp):
    #parallelizable function to tabix final dose and empiricalDose vcfs

    result1 = subprocess.call(f'tabix -f -p vcf {fp}', shell=True)
    return result1


def write_vcf_pos_to_file(vcf, outfp):
    #parse vcf file for positions and write to outfp
    result1 = subprocess.call(f"zcat {vcf} | grep ^[^#] | cut -f2 > {outfp}", shell=True)


#
# Start
#

if __name__ == '__main__':

    #Init. log and folders
    write_log(log,f'target vcf: {tar}\n')
    append_log(log,f'nophasing set to {nophasing}\n')


    Path(outdir).mkdir(parents=True, exist_ok=True)
    Path(qcdir).mkdir(parents=True, exist_ok=True)
    Path(chr_dir).mkdir(parents=True, exist_ok=True)
    Path(impute_chr_dir).mkdir(parents=True, exist_ok=True)


    #Extract target file variant positions
    append_log(log,f'creating and parsing pos file: {posfp}\n')
    write_vcf_pos_to_file(tar, posfp)
    with open(posfp) as f:
        positions = [int(line) for line in f.read().splitlines()]


    #Make QC folders
    append_log(log,f'creating output directories\n')
    chrfilter_dir = f'{qcdir}/{chrom}_filteredChunk'
    Path(chr_dir).mkdir(parents=True, exist_ok=True)


    #Create windows
    #Take first position and get closest 1Mb + 1 (e.g. 1, 1000001, 2000001, 3000001, etc.)
    append_log(log,f'creating windows:\n')
    startfloor = positions[0]//min_window_offset
    startpos = startfloor*min_window_offset + 1
    windows = position_windows(pos=np.array(positions), size=int(chunksize), start=startpos, step=int(chunksize-overlap))
    append_log(log,f'{windows}\n')


    #Create imputedir
    append_log(log,f'reading positions from ref: {ref}\n')
    ref_pd, ref_pos = get_chrompos_pd_from_vcf(ref)
    with open(f'{outdir}/{chrom}.ref.pos.txt', 'w') as g:
        for pos in ref_pos:
            g.write(f'{pos}\n')


    #Extract chunks
    append_log(log,f'creating chunks\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        create_chunks_results = executor.map(create_chunks, windows)
        executor.shutdown()
    for result, window in zip(create_chunks_results, windows):
        append_log(log,f'create_chunks on {window}: {result}\n')
    chunk_glob = glob.glob(chr_dir+f'/chunk_*_*.vcf.gz')
    append_log(log,f'expected: {len(windows)}, found: {len(chunk_glob)}\n')


    #QC chunks
    append_log(log,f'qc chunks\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        qc_results = executor.map(qc_chunks, windows)
        executor.shutdown()
    for result, window in zip(qc_results, windows):
        append_log(log,f'qc on {window}: {result}\n')
    qc_glob = glob.glob(f'{chrfilter_dir}/chunk_*_*.txt')
    append_log(log,f'expected: {len(windows)}, found: {len(qc_glob)}\n')


    #Update windows based on QC results
    append_log(log,f'checking windows after QC\n')
    prev_windows = list(windows)
    windows = []
    for window in prev_windows:
        w_dict = create_pathdict_by_window(window)
        filteredchunk = open(w_dict['chrfilterchunk']).read().strip()
        if filteredchunk == 'Keep':
            windows.append(window)
        else:
            append_log(log,f'removing window {window} from analysis: {filteredchunk}\n')


    #BCF/filtering chunks
    append_log(log,f'bcf chunks\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        filter_results = executor.map(filter_chunks, windows)
        executor.shutdown()
    for result, window in zip(filter_results, windows):
        append_log(log,f'filter on {window}: {result}\n')

    bcf_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.bcf')
    append_log(log,f'expected: {len(windows)}, found: {len(bcf_glob)}\n')


    #Phasing chunks unless opt out
    if nophasing:
        #convert bcf to vcf for imputation step
        append_log(log,f'converting chunks (without phasing)\n')
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
            vcf_results = executor.map(convert_bcf_to_vcf_without_phasing, windows)
            executor.shutdown()
        phase_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.notrephased.vcf.gz')
        append_log(log,f'expected: {len(windows)}, found: {len(phase_glob)}\n')
    else:
        #with multicore workers each
        append_log(log,f'phase chunks\n')
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
            vcf_results = executor.map(phase_chunks, windows)
            executor.shutdown()
        phase_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.phased.vcf.gz')
        append_log(log,f'expected: {len(windows)}, found: {len(phase_glob)}\n')
    for result, window in zip(vcf_results, windows):
        append_log(log,f'vcf on {window}: {result}\n')


    #Create imputation dir for chrom
    Path(impute_chr_dir).mkdir(parents=True, exist_ok=True)


    #Impute chunks with multicore workers
    append_log(log,f'impute chunks\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
        impute_results = executor.map(impute_chunks, windows)
        executor.shutdown()
    for result, window in zip(impute_results, windows):
        append_log(log,f'impute on {window}: {result}\n')
    impute_glob = glob.glob(f'{impute_chr_dir}/chunk_*_*.imputed.dose.vcf.gz')
    append_log(log,f'expected: {len(windows)}, found: {len(impute_glob)}\n')


    #Gather imputed chunk filepaths
    topmed_imputation_aggregate_chunks = [f'{impute_chr_dir}/chunk_{start}_{stop}.imputed.dose.vcf.gz' for start,stop in windows]
    for chunk in topmed_imputation_aggregate_chunks:
        assert os.path.exists(f'{chunk}') == True, '%s does not exist' % (chunk)
    append_log(log,'topmed_imputation_aggregate_chunks:\n')
    for chunk in topmed_imputation_aggregate_chunks:
        append_log(log,chunk+'\n')


    #Write summary for chunks
    filtered_chunk_vcfs = []
    chr_summary_file = f'{outdir}/{chrom}_chunk_summary.txt'
    append_log(log,f'Write summary of each chunk to {chr_summary_file}\n')
    with open(chr_summary_file, 'w') as g:
        for f in topmed_imputation_aggregate_chunks:
            append_log(log,f'parsing input: {f} ({os.stat(f).st_size})\n')
            if os.stat(f).st_size != 0:
                start = int(os.path.basename(f).split('_')[1])
                filtered_chunk_vcfs.append((start, f))
                g.write(f'{f}\tKeep\n')
            else:
                g.write(f'{f}\tSkip (QC or no valid variants to impute)\n')


    #Get final trim windows from vcfs
    append_log(log,'Get final chunk windows from vcf path:\n')
    chunk_vcfs = pd.DataFrame(filtered_chunk_vcfs, columns=['start','file'])
    chunk_vcfs = chunk_vcfs.sort_values('start')
    chunk_vcfs = list(chunk_vcfs['file'])
    final_chunk_windows = process_chunk_windows_to_final_trim_windows(chunk_vcfs)
    for fp,final_window in zip(chunk_vcfs,final_chunk_windows):
        append_log(log,f'{fp}: {final_window}\n')


    #Final output file names/paths
    imputchromvcf = f'{outdir}/{chrom}.dose.vcf.gz'
    imputchromvcfED = f'{outdir}/{chrom}.empiricalDose.vcf.gz'
    imputchrominfo = f'{outdir}/{chrom}.info.gz'


    #Concat final dose file
    append_log(log,'Generating final dose vcf\n')
    final_dose_result = write_final_vcf_file(chunk_vcfs, final_chunk_windows, imputchromvcf)


    #Concat final empiricalDose file
    append_log(log,'Generating final empiricalDose vcf\n')
    chunk_vcfs_ed = [chunk.replace('.dose.vcf.gz', '.empiricalDose.vcf.gz') for chunk in chunk_vcfs]
    final_empericaldose_result = write_final_vcf_file(chunk_vcfs_ed, final_chunk_windows, imputchromvcfED)


    #Write final info file
    append_log(log,'Writing final info\n')
    chunk_infos = [chunk.replace('.dose.vcf.gz', '.info') for chunk in chunk_vcfs]
    write_final_info_file(chunk_infos, final_chunk_windows, imputchrominfo)


    #Clean up temp. files on success
    if final_dose_result == 0 and final_empericaldose_result == 0:
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)


    #Tabix dose and empiricalDose vcfs in parallel
    append_log(log,'Indexing final dose and empiricalDose vcfs\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        tabix_final_results = executor.map(HELPER_tabix_vcf_parallel, [imputchromvcf, imputchromvcfED])
        executor.shutdown()
    output = subprocess.check_output(f'bcftools index -n {imputchromvcf}', shell=True)
    append_log(log,f'number of snps in dosage file: {output.decode("utf-8")}\n')


    #Make final histogram: # imputed variants vs. # of variants in panel
    if not skiphistplot:
        tar_pos_fp = f'{outdir}/{chrom}.tar.pos.txt'
        append_log(log,f'Writing target pos to {tar_pos_fp}\n')
        write_vcf_pos_to_file(imputchromvcf, tar_pos_fp)
        append_log(log,f'Plotting target vs ref histogram\n')
        subprocess.call(f'Rscript {HIST_R_SCRIPT} --ref {outdir}/{chrom}.ref.pos.txt \
                                                  --tar {outdir}/{chrom}.tar.pos.txt \
                                                  --out {outdir}/{chrom}.tar.vs.ref.hist.pdf \
                                                  --prefix "{chrom} tar vs. ref"', shell=True)


