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

# based on implementation of the TOPMed pipeline by Minhui Chen

#
# Program paths
#

# This script will also use bcftools and tabix!

PLINK = '/project/haiman_625/Software/imputation_pipeline_CharlestonGroup/programs/plink'
EAGLE = '/project/haiman_625/Software/imputation_pipeline_CharlestonGroup/programs/eagle'
GMAP = '/project/haiman_625/Software/imputation_pipeline_CharlestonGroup/programs/genetic_map_hg38_withX.txt.gz'
MM4 = '/project/chia657_28/programs/minimac4/bin/minimac4'

parser = argparse.ArgumentParser()
parser.add_argument('--tar', required=True, help='target data VCF')
parser.add_argument('--ref', required=True, help='phased reference VCF')
parser.add_argument('--refm3', required=True, help='reference m3vcf')
parser.add_argument('--chrom', required=True, help='format the same way as hg38 (e.g. "chr22")')
parser.add_argument('--outdir', required=True, help='without trailing /')
parser.add_argument('--workers', required=True, help='number of chunks to work on at same time')
parser.add_argument('--mem', required=True, help='mem in gb to use')
parser.add_argument('--windowsize', required=False, help='mm4 window size in bp', default=500000)
parser.add_argument('--chunksize', required=False, help='chunk size in bp', default=25e6)
parser.add_argument('--overlap', required=False, help='overlap size in bp', default=5e6)
parser.add_argument('--multicore', required=False, help='number of workers to use when multithreading', default=5)
parser.add_argument('--nophasing', default=False, action='store_true', help='Skips phasing during pipeline (for analyses such as meta-imputation)')

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

print(nophasing)

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
buffersize = 2.5e6

#
# Functions
#


def position_windows(pos, size, start=None, stop=None, step=None):
    """Convenience function to construct windows for the
    :func:`windowed_statistic` and :func:`windowed_count` functions.
    Params:
    -------
    pos: array-like, not series
    """
    last = False
    # determine start and stop positions
    if start is None:
        start = pos[0]
    if stop is None:
        stop = pos[-1]
    if step is None:
        # non-overlapping
        step = size
    windows = []
    for window_start in range(start, stop, step):
        # determine window stop
        window_stop = window_start + size
        if window_stop >= stop:
            # last window
            last = True
        else:
            window_stop -= 1
        windows.append([window_start, window_stop])
        if last:
            break
    return np.asarray(windows)


def get_chrompos_pd_from_vcf(filepath):
    #read in vcf and parse for chrom, pos, ref, and alt

    if filepath.endswith('.gz'):
        cat_prefix = 'z'
    else:
        cat_prefix = ''
    lines = [l.split() for l in subprocess.check_output(f"{cat_prefix}cat {filepath} | grep ^[^#] | cut -f1-2,4-5", shell=True).decode().splitlines()]
    lines = [f'{l[0]},{l[1]},{sorted([l[2],l[3]])[0]}{sorted([l[2],l[3]])[1]}' for l in lines]
    df = pd.DataFrame(lines, columns=['snp'])
    return df


def create_chunks(window):
    #parse vcf between start and end
    start, stop = window
    chunkprefix = chr_dir+f'/chunk_{start}_{stop}'

    print('starting create_chunks',window)
    #subprocess.call(f'bcftools view -Oz -o {chunkprefix}.vcf.gz -r {chrom}:{start}-{stop} {tar}', shell=True, stdout=subprocess.DEVNULL)
    subprocess.call(f'bcftools view -Oz -o {chunkprefix}.vcf.gz -r {chrom}:{start}-{stop} {tar}', shell=True)
    print('finished create_chunks',window)


def create_pathdict_by_window(window):
    start, stop = window

    qcdir = outdir+'/qc'
    chr_dir = f'{qcdir}/{chrom}'
    impute_dir = outdir+f'/imputation/{chrom}'
    d = {}

    d['chunkprefix'] = chr_dir+f'/chunk_{start}_{stop}'
    d['chunkvcf'] = chr_dir+f'/chunk_{start}_{stop}.vcf.gz'
    d['chrmissing_dir'] = f'{qcdir}/{chrom}_missing'

    d['chunkmissingprefix'] = f'{d["chrmissing_dir"]}/chunk_{start}_{stop}'
    d['chrfilter_dir'] = f'{qcdir}/{chrom}_filteredChunk'
    d['chrfilterchunk'] = f'{d["chrfilter_dir"]}/chunk_{start}_{stop}.txt'

    #d['filteredprefix'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}'
    #d['filteredprefixwithchr'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr'
    d['snpsfile'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.snps'

    d['bcf'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.bcf'
    d['index'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.bcf.csi'

    d['chrmissing_dir'] = f'{qcdir}/{chrom}_missing'
    d['chunkmissingprefix'] = f'{d["chrmissing_dir"]}/chunk_{start}_{stop}'

    d['phasedprefix'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.phased'
    d['phasedvcf'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.phased.vcf.gz'
    d['skipphasevcf'] = d['chrfilter_dir']+f'/chunk_{start}_{stop}.withchr.notrephased.vcf.gz'

    d['imputeddosefile'] = impute_dir+f'/chunk_{start}_{stop}.imputed.dose.vcf.gz'
    d['imputedprefix'] = impute_dir+f'/chunk_{start}_{stop}.imputed'

    return d


def qc_chunks(window):
    print('starting qc_chunks',window)
    w_dict = create_pathdict_by_window(window)

    chunk_pd = get_chrompos_pd_from_vcf(w_dict['chunkvcf'])
    num_valid_variants = ref_pd.loc[ref_pd['snp'].isin(chunk_pd['snp'])].shape[0]
    percent_valid_variants = num_valid_variants/ref_pd.shape[0]
    Path(w_dict['chrmissing_dir']).mkdir(parents=True, exist_ok=True)
    subprocess.call(f'{PLINK} --vcf {w_dict["chunkvcf"]} --missing --memory {mem_per_worker} --out {w_dict["chunkmissingprefix"]}', shell=True, stdout=subprocess.DEVNULL)
    missing = pd.read_csv(w_dict['chunkmissingprefix']+'.imiss', sep='\s+')
    num_high_missing = missing.loc[missing['F_MISS'] > (1-minsampleCallRate)].shape[0]
    Path(w_dict['chrfilter_dir']).mkdir(parents=True, exist_ok=True)

    with open(w_dict['chrfilterchunk'], 'w') as f:
        if num_valid_variants < min_valid_variants or num_high_missing > 0:
            f.write(f'Exclude\t{num_valid_variants}\t{percent_valid_variants}\t{num_high_missing}')
        else:
            f.write('Keep')
    print('finished qc_chunks',window)



def filter_chunks(window):
    w_dict = create_pathdict_by_window(window)

    filteredchunk = open(w_dict['chrfilterchunk']).read().strip()

    if filteredchunk == 'Keep':
        print(f'filter_chunks == keep:',window)
        #read bim into df
        lines = [l.split() for l in subprocess.check_output(f"zcat {w_dict['chunkvcf']} | grep ^[^#] | cut -f1-5", shell=True).decode().splitlines()]
        print(f'finished read chunkvcf into list',window)
        lines = [[l[0], l[1], l[2], sorted([l[3],l[4]])[0], sorted([l[3],l[4]])[1]] for l in lines]
        bim = pd.DataFrame(lines, columns=['chr', 'pos', 'origsnpid','a1', 'a2'])
        print(f'finished create bim df from list',window)
        bim['pos'] = bim['pos'].astype('Int64')
        bim['snp'] = bim.apply(lambda row: f"{row['chr']},{row['pos']},{row['a1']}{row['a2']}", axis=1)
        print(f'finished create chrompos col in df',window)

        # pick loc of snps only
        bim = bim.loc[(bim['a1'].isin(alleles+['0'])) & (bim['a2'].isin(alleles+['0']))]
        pos = np.array(bim['pos'])
        print(f'finished subset pos in df',window)

        # make sure dis is not negative
        dis = pos[1:] - pos[:-1]
        if (dis < 0).sum() > 0:
            print(f'Positions are not inc in bim')
            return

        snp = np.array(bim['snp'])
        duplicates = np.append(snp[:-1][dis<=1], snp[1:][dis<=1])
        bim = bim.loc[~bim['snp'].isin(duplicates)]

        print(f'finished remove dup pos in df',window)

        # filter missing rate / call rate
        lmissfile = w_dict['chunkmissingprefix']+'.lmiss'
        highmissing = pd.read_csv(lmissfile, sep='\s+')
        print(f'finished parsing lmissfile',window)
        highmissing = highmissing.loc[highmissing['F_MISS'] > (1-minsnpCallRate)]
        bim = bim.loc[~bim['snp'].isin(highmissing['SNP'])]
        print(f'finished remove highmissing in df',window)

        # extract snps and convert to bcf, then index
        #bim['snp'].to_csv(w_dict['snpsfile'], sep='\t', index=False, header=False)
        bim['origsnpid'].to_csv(w_dict['snpsfile'], sep='\t', index=False, header=False)
        print(w_dict['snpsfile'])
        bcftools_cmd = f'bcftools view --include ID==@{w_dict["snpsfile"]} -Ob -o {w_dict["bcf"]} {w_dict["chunkvcf"]}'
        subprocess.call(bcftools_cmd, shell=True)
        print(f'finished bcftools call',window)
        subprocess.call(f'bcftools index {w_dict["bcf"]}', shell=True)
        print(f'finished bcftools indexing',window)
    else:
        print(f'filter_chunks != keep: touch {w_dict["bcf"]}',window)
        subprocess.call(f'touch {w_dict["bcf"]} {w_dict["index"]}', shell=True)


def phase_chunks(window):
    w_dict = create_pathdict_by_window(window)
    start, stop = window

    filteredchunk = open(w_dict['chrfilterchunk']).read().strip()
    if filteredchunk == 'Keep':
        subprocess.call(f'{EAGLE} --vcfRef {ref} \
                                  --vcfTarget {w_dict["bcf"]} \
                                  --geneticMapFile {GMAP} \
                                  --chrom {chrom} \
                                  --bpStart {start} \
                                  --bpEnd {stop} \
                                  --outPrefix {w_dict["phasedprefix"]} \
                                  --numThreads {multicore} \
                                  --allowRefAltSwap \
                                  --vcfOutFormat z', shell=True, stdout=subprocess.DEVNULL)
    else:
        subprocess.call(f'touch {w_dict["phasedvcf"]}', shell=True)


def convert_bcf_to_vcf_without_phasing(window):
    # previously plink -> bcf -> vcf
    # now vcf -> bcf -> vcf

    w_dict = create_pathdict_by_window(window)

    filteredchunk = open(w_dict['chrfilterchunk']).read().strip()
    if filteredchunk == 'Keep':
        subprocess.call(f'bcftools view -Oz -o {w_dict["skipphasevcf"]} --threads {multicore} {w_dict["bcf"]}', shell=True, stdout=subprocess.DEVNULL)
    else:
        subprocess.call(f'touch {w_dict["skipphasevcf"]}', shell=True)


def impute_chunks(window):
    w_dict = create_pathdict_by_window(window)
    start, stop = window

    if nophasing:
        preimputationvcf = w_dict['skipphasevcf']
    else:
        preimputationvcf = w_dict['phasedvcf']

    filteredchunk = open(w_dict['chrfilterchunk']).read().strip()
    print(filteredchunk)

    if filteredchunk == 'Keep' and stop - buffersize > start + buffersize:
        # check whether variants exist in imputation region:
        # chunk 1:120000001-145000000 pass QC filter, but don't have variants in
        # imputation region 1:122000001-143000000
        empty_imputation_region = True
        for line in gzip.open(preimputationvcf, 'rt'):
            if not re.search('^#', line):
                line = line.split('\t', 2)
                pos = int(line[1])
                if pos > start + buffersize and pos < stop - buffersize:
                    empty_imputation_region = False
                    break
        if empty_imputation_region:
            cmd = f'touch {w_dict["imputeddosefile"]}'
            print(f'empty impute region, running:{cmd}')
            subprocess.call(cmd, shell=True)
        else:
            subprocess.call(f'{MM4} --refHaps {refm3} \
                                    --haps {preimputationvcf} --chr {chrom} \
                                    --start {start} \
                                    --end {stop} \
                                    --window {windowsize} \
                                    --cpus {multicore} \
                                    --format GT,DS,GP \
                                    --allTypedSites \
                                    --meta \
                                    --noPhoneHome \
                                    --minRatio 0.00001 \
                                    --prefix {w_dict["imputedprefix"]}', shell=True)
            if os.path.exists(f'{w_dict["imputeddosefile"]}.tbi'):
                os.remove(f'{w_dict["imputeddosefile"]}.tbi')
            subprocess.call(f'tabix -f -p vcf {w_dict["imputeddosefile"]}', shell=True, stdout=subprocess.DEVNULL)
            if os.path.exists(f'{w_dict["imputeddosefile"].replace("dose","empiricalDose")}.tbi'):
                os.remove(f'{w_dict["imputeddosefile"].replace("dose","empiricalDose")}.tbi')
            subprocess.call(f'tabix -f -p vcf {w_dict["imputeddosefile"].replace("dose","empiricalDose")}', shell=True, stdout=subprocess.DEVNULL)
    else:
        cmd = f'touch {w_dict["imputeddosefile"]}'
        print(f'failed "Keep" and start/stop check, running:{cmd}')
        subprocess.call(cmd, shell=True)


def parse_info_between_pos(info_path, start_bp, end_bp):
    #read in an info file and only keep the lines between the interval (inclusive) based on SNP name

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


#
# Start
#

if __name__ == '__main__':

    # create outfolder
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # calculate intervals for chunks
    posdir = outdir+'/targetpositions'
    Path(posdir).mkdir(parents=True, exist_ok=True)
    pos = posdir+f'/{chrom}.pos.txt'
    temp_pos = pos+'.temp'

    #write to temp file and move if completed/successful
    cat_prefix = 'z' if tar.endswith('.gz') else ''
    subprocess.run(f"{cat_prefix}cat {tar} | grep ^[^#] | cut -f2 > {temp_pos}", shell=True)
    subprocess.run(f"mv {temp_pos} {pos}", shell=True)

    with open(pos) as f:
        positions = [int(line) for line in f.read().splitlines()]

    #qc target data
    qcdir = outdir+'/qc'
    Path(qcdir).mkdir(parents=True, exist_ok=True)

    #extract target chunks
    chr_dir = f'{qcdir}/{chrom}'
    chrfilter_dir = f'{qcdir}/{chrom}_filteredChunk'
    imputdir = outdir+'/imputation'
    if os.path.isdir(chr_dir):
        shutil.rmtree(chr_dir)
    Path(chr_dir).mkdir(parents=True, exist_ok=True)

    windows = position_windows(pos=np.array(positions), size=int(chunksize), start=1, step=int(chunksize-overlap))

    Path(imputdir).mkdir(parents=True, exist_ok=True)

    ref_pd = get_chrompos_pd_from_vcf(ref)

    #with 1 worker each
    print(f'creating chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(create_chunks, windows)
        executor.shutdown()

    chunk_glob = glob.glob(chr_dir+f'/chunk_*_*.vcf.gz')
    print(f'expected: {len(windows)}, found: {len(chunk_glob)}')



    print(f'qc chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(qc_chunks, windows)
        executor.shutdown()


    qc_glob = glob.glob(f'{chrfilter_dir}/chunk_*_*.txt')
    print(f'expected: {len(windows)}, found: {len(qc_glob)}')



    print(f'bcf chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(filter_chunks, windows)
        executor.shutdown()
    bcf_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.bcf')
    print(f'expected: {len(windows)}, found: {len(bcf_glob)}')



    if nophasing:
        #convert bcf to vcf for imputation step
        print(f'converting chunks (without phasing)')
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
            executor.map(convert_bcf_to_vcf_without_phasing, windows)
            executor.shutdown()
        phase_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.notrephased.vcf.gz')
        print(f'expected: {len(windows)}, found: {len(phase_glob)}')
    else:
        #with multicore workers each
        print(f'phase chunks')
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
            executor.map(phase_chunks, windows)
            executor.shutdown()
        phase_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.phased.vcf.gz')
        print(f'expected: {len(windows)}, found: {len(phase_glob)}')


    impute_dir = outdir+f'/imputation/{chrom}'
    Path(impute_dir).mkdir(parents=True, exist_ok=True)

    #with multicore workers each
    print(f'impute chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
        executor.map(impute_chunks, windows)
        executor.shutdown()
    impute_glob = glob.glob(f'{outdir}/imputation/{chrom}/chunk_*_*.imputed.dose.vcf.gz')
    print(f'expected: {len(windows)}, found: {len(impute_glob)}')

    topmed_imputation_aggregate_chunks = [f'{outdir}/imputation/{chrom}/chunk_{start}_{stop}.imputed.dose.vcf.gz' for start,stop in windows]
    for chunk in topmed_imputation_aggregate_chunks:
        #assert len(windows) == len(impute_glob), 'num imputed chunks (%s) != num windows (%s) ' % (len(windows), len(impute_glob))
        assert os.path.exists(f'{chunk}') == True, '%s does not exist' % (chunk)



    #stitch imputed chunks together
    outchunks = f'{imputdir}/{chrom}.imputed.dose.vcf.chunks'

    #dtype = [('start', int), ('file', 'U100')]
    chunk_vcfs = []

    print('topmed_imputation_aggregate_chunks:')
    print(topmed_imputation_aggregate_chunks)

    #write summary file of each window, e.g.:
    chr_summary_file = f'{imputdir}/{chrom}_chunk_summary.txt'
    with open(chr_summary_file, 'w') as g:

        for f in topmed_imputation_aggregate_chunks:
            print(f'parsing input: {f} ({os.stat(f).st_size})')
            if os.stat(f).st_size != 0:
                start = int(os.path.basename(f).split('_')[1])
                chunk_vcfs.append((start, f))
                g.write(f'chunk_{start}_{stop}\tKeep\n')
            else:
                g.write(f'chunk_{start}_{stop}\tSkip\n')


    chunk_vcfs = pd.DataFrame(chunk_vcfs, columns=['start','file'])
    chunk_vcfs = chunk_vcfs.sort_values('start')
    chunk_vcfs = list(chunk_vcfs['file'])

    imputchromvcf = f'{imputdir}/{chrom}.dose.vcf.gz'
    imputchromvcfED = f'{imputdir}/{chrom}.empiricalDose.vcf.gz'
    imputchrominfo = f'{imputdir}/{chrom}.info.gz'

    new_chunk_vcfs = []

    #check if next file overlaps
    idx = 0
    if len(chunk_vcfs) > 0:
        chunk_infos = [chunk.replace('dose.vcf.gz','info') for chunk in chunk_vcfs]

        with gzip.open(imputchrominfo, 'wt') as g:
            #write info header
            g.write('SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1\n')
            while idx < len(chunk_vcfs)-1:
                chunk = chunk_vcfs[idx]
                chunkED = chunk.replace('dose','empiricalDose')
                info = chunk_infos[idx]
                next_chunk = chunk_vcfs[idx+1]
                next_chunkED = next_chunk.replace('dose','empiricalDose')

                if os.path.exists(f'{chunk}.tbi'):
                    os.remove(f'{chunk}.tbi')

                if os.path.exists(f'{chunkED}.tbi'):
                    os.remove(f'{chunkED}.tbi')

                print(f'tabix -f -p vcf {chunk}')
                subprocess.run(f'tabix -f -p vcf {chunk}', shell=True)
                print(f'tabix -f -p vcf {chunkED}')
                subprocess.run(f'tabix -f -p vcf {chunkED}', shell=True)

                if os.path.exists(f'{next_chunk}.tbi'):
                    os.remove(f'{next_chunk}.tbi')
                if os.path.exists(f"{next_chunkED}.tbi"):
                    os.remove(f"{next_chunkED}.tbi")

                print(f'tabix -f -p vcf {next_chunk}')
                subprocess.run(f'tabix -f -p vcf {next_chunk}', shell=True)
                print(f"tabix -f -p vcf {next_chunkED}")
                subprocess.run(f"tabix -f -p vcf {next_chunkED}", shell=True)

                print(chunk)

                cur_start = int(os.path.basename(chunk).split('_')[1])
                cur_stop = int(os.path.basename(chunk).split('_')[2].split('.')[0])
                next_start = int(os.path.basename(next_chunk).split('_')[1])
                next_stop = int(os.path.basename(next_chunk).split('_')[2].split('.')[0])

                if next_start < cur_stop:
                    #overlap!

                    #change cur stop to new pos and extract vcf, update to new path
                    new_cur_stop = cur_stop - int(overlap/2)
                    new_chunk = chunk.replace(str(cur_stop), str(new_cur_stop))
                    new_chunk_vcfs.append(new_chunk)
                    new_chunkED = new_chunk.replace('dose','empiricalDose')

                    subprocess.call(f"bcftools view -r {chrom}:{cur_start}-{new_cur_stop} -Oz -o {new_chunk} {chunk}", shell=True)
                    if os.path.isfile(new_chunk) and os.path.getsize(new_chunk) > 0:
                        os.remove(chunk)

                    subprocess.call(f"bcftools view -r {chrom}:{cur_start}-{new_cur_stop} -Oz -o {new_chunkED} {chunkED}", shell=True)
                    if os.path.isfile(new_chunkED) and os.path.getsize(new_chunkED) > 0:
                        os.remove(chunkED)

                    if os.path.exists(f'{new_chunk}.tbi'):
                        os.remove(f'{new_chunk}.tbi')
                    if os.path.exists(f'{new_chunkED}.tbi'):
                        os.remove(f'{new_chunkED}.tbi')

                    print(f'tabix -f -p vcf {new_chunk}')
                    subprocess.check_output(f'tabix -f -p vcf {new_chunk}', shell=True)
                    print(f'tabix -f -p vcf {new_chunkED}')
                    subprocess.check_output(f'tabix -f -p vcf {new_chunkED}', shell=True)
                    #chunk_vcfs[idx] = new_chunk

                    #write info up to new stop
                    valid_info_lines = parse_info_between_pos(info, cur_start, new_cur_stop)
                    for line in valid_info_lines:
                        g.write(f'{line}\n')

                    #change next_start to new pos, extract vcf, update to new path
                    new_next_start = next_start + int(overlap/2)
                    new_next_chunk = next_chunk.replace(str(next_start), str(new_next_start))
                    new_next_chunkED = new_next_chunk.replace('dose','empiricalDose')
                    subprocess.call(f'bcftools view -r {chrom}:{new_next_start}-{next_stop} -Oz -o {new_next_chunk} {next_chunk}', shell=True)
                    if os.path.isfile(new_next_chunk) and os.path.getsize(new_next_chunk) > 0:
                        os.remove(next_chunk)
                    subprocess.call(f'bcftools view -r {chrom}:{new_next_start}-{next_stop} -Oz -o {new_next_chunkED} {next_chunkED}', shell=True)
                    if os.path.isfile(new_next_chunkED) and os.path.getsize(new_next_chunkED) > 0:
                        os.remove(next_chunkED)

                    if os.path.exists(f'{new_next_chunk}.tbi'):
                        os.remove(f'{new_next_chunk}.tbi')
                    if os.path.exists(f'{new_next_chunkED}.tbi'):
                        os.remove(f'{new_next_chunkED}.tbi')
                    print(f'tabix -f -p vcf {new_next_chunk}')
                    subprocess.check_output(f'tabix -f -p vcf {new_next_chunk}', shell=True)
                    print(f'tabix -f -p vcf {new_next_chunkED}')
                    subprocess.check_output(f'tabix -f -p vcf {new_next_chunkED}', shell=True)
                    chunk_vcfs[idx+1] = new_next_chunk
                else:
                    #no overlap
                    #can write full info for current chunk
                    with open(info) as f:
                        lines = [rline.strip() for rline in f.read().splitlines()[1:]]
                    for line in lines:
                        g.write(f'{line}\n')
                idx += 1
            if idx == len(chunk_vcfs)-1:
                chunk = chunk_vcfs[idx]
                new_chunk_vcfs.append(chunk)
                info = chunk_infos[idx]
                cur_start = int(os.path.basename(chunk).split('_')[1])
                cur_stop = int(os.path.basename(chunk).split('_')[2].split('.')[0])

                valid_info_lines = parse_info_between_pos(info, cur_start, cur_stop)
                for line in valid_info_lines:
                    g.write(f'{line}\n')

    print(f'writing to {outchunks}')
    with open(outchunks, 'w') as f:
        f.write('\n'.join(new_chunk_vcfs))
    print(f'writing done')

    chunkstring = ' '.join(new_chunk_vcfs)
    chunkstringED = ' '.join([chunk.replace('dose','empiricalDose') for chunk in new_chunk_vcfs])
    print('concatenating chunks')

    subprocess.call(f'bcftools concat --threads {workers} -Oz -o {imputchromvcf} -a -d all {chunkstring}', shell=True)
    if os.path.isfile(imputchromvcf) and os.path.getsize(imputchromvcf) > 0:
        for fp in chunkstring.split():
            os.remove(fp)

    subprocess.call(f'bcftools concat --threads {workers} -Oz -o {imputchromvcfED} -a -d all {chunkstringED}', shell=True)
    if os.path.isfile(imputchromvcfED) and os.path.getsize(imputchromvcfED) > 0:
        for fp in chunkstringED.split():
            os.remove(fp)

    if os.path.exists(f'{imputchromvcf}.tbi'):
        os.remove(f'{imputchromvcf}.tbi')
    subprocess.run(f'tabix -f -p vcf {imputchromvcf}', shell=True)
    if os.path.exists(f'{imputchromvcfED}.tbi'):
        os.remove(f'{imputchromvcfED}.tbi')
    subprocess.run(f'tabix -f -p vcf {imputchromvcfED}', shell=True)
    print('concatenating and tabix done')


