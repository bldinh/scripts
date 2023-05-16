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

PLINK = '/scratch1/bldinh/programs/plink'
EAGLE = '/scratch1/bldinh/programs/eagle'
GMAP = '/scratch1/bldinh/programs/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz'
MM4 = '/project/chia657_28/programs/minimac4/bin/minimac4'

parser = argparse.ArgumentParser()
parser.add_argument('--tar', required=True, help='prefix for plink format target data')
parser.add_argument('--ref', required=True, help='phased reference data')
parser.add_argument('--refm3', required=True, help='m3vcf')
parser.add_argument('--refbim', required=True, help='bim')
parser.add_argument('--chrom', required=True, help='format the same way as hg38 (e.g. "chr22")')
parser.add_argument('--outdir', required=True, help='without trailing /')
parser.add_argument('--workers', required=True, help='number of chunks to work on at same time')
parser.add_argument('--mem', required=True, help='mem in gb to use')
parser.add_argument('--windowsize', required=False, help='mm4 window size in bp', default=500000)
parser.add_argument('--chunksize', required=False, help='chunk size in bp', default=25e6)
parser.add_argument('--overlap', required=False, help='overlap size in bp', default=5e6)
parser.add_argument('--multicore', required=False, help='number of workers to use when multithreading', default=5)

args = parser.parse_args()
tar = args.tar
ref = args.ref
refm3 = args.refm3
refbim = args.refbim
chrom = args.chrom
outdir = args.outdir
workers = int(args.workers)
mem = float(args.mem)
chunksize = int(args.chunksize)
overlap = int(args.overlap)
windowsize = int(args.windowsize)
multicore = int(args.multicore)

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
            window_stop = stop
            last = True
        else:
            window_stop -= 1
        windows.append([window_start, window_stop])
        if last:
            break
    return np.asarray(windows)


def update_bim_snpname(bim_fn):
    '''
    Update bim file snpname to chr_pos_a2_a1 (a2 a1 are ref alt if keep-allele-order from vcf)
    Parameters
    ----------
    bim_fn : str
        filename of bim
    Returns
    -------
    Notes
    -----
    Replace the original bim file.
    '''
    bim = pd.read_csv(bim_fn, sep='\t', header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2']) # a2 a1 are ref alt
    bim['snp'] = bim['chr'].astype('str')+'_'+bim['pos'].astype('str')+'_'+bim['a2']+'_'+bim['a1']
    bim.to_csv(bim_fn, sep='\t', index=False, header=False)


def create_chunks(window):
    #parse vcf between start and end

    start, stop = window
    chunkprefix = chr_dir+f'/chunk_{start}_{stop}'

    subprocess.call(f'{PLINK} --bfile {tar} --chr {chrom} --from-bp {start} --to-bp {stop} --memory {mem_per_worker} --make-bed --out {chunkprefix}', shell=True, stdout=subprocess.DEVNULL)
    update_bim_snpname(chunkprefix+'.bim')


def qc_chunks(window):
    #parse vcf between start and end

    start, stop = window
    chunkprefix = chr_dir+f'/chunk_{start}_{stop}'
    chunkbim = chr_dir+f'/chunk_{start}_{stop}.bim'
    chunkbi_pd = pd.read_csv(chunkbim, sep='\s+', header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2'])

    no_valid_variants = refbim_pd.loc[refbim_pd['snp'].isin(chunkbi_pd['snp'])].shape[0]
    percent_valid_variants = no_valid_variants/refbim_pd.shape[0]

    chrmissing_dir = f'{qcdir}/{chrom}_missing'
    Path(chrmissing_dir).mkdir(parents=True, exist_ok=True)

    chunkmissingprefix = f'{chrmissing_dir}/chunk_{start}_{stop}'
    subprocess.call(f'{PLINK} --bfile {chunkprefix} --missing --memory {mem_per_worker} --out {chunkmissingprefix}', shell=True, stdout=subprocess.DEVNULL)
    missing = pd.read_csv(chunkmissingprefix+'.imiss', sep='\s+')
    no_high_missing = missing.loc[missing['F_MISS'] > (1-minsampleCallRate)].shape[0]

    chrfilter_dir = f'{qcdir}/{chrom}_filteredChunk'
    Path(chrfilter_dir).mkdir(parents=True, exist_ok=True)
    chrfilterchunk = f'{chrfilter_dir}/chunk_{start}_{stop}.txt'

    with open(chrfilterchunk, 'w') as f:
        #if no_valid_variants < min_valid_variants or percent_valid_variants < min_percent_valid_variants or no_high_missing > 0:
        if no_valid_variants < min_valid_variants or no_high_missing > 0:
            f.write(f'Exclude\t{no_valid_variants}\t{percent_valid_variants}\t{no_high_missing}')
        else:
            f.write('Keep')


def filter_chunks(window):
    start, stop = window
    chunkprefix = chr_dir+f'/chunk_{start}_{stop}'
    chunkbim = chr_dir+f'/chunk_{start}_{stop}.bim'
    chrfilter_dir = f'{qcdir}/{chrom}_filteredChunk'
    chrfilterchunk = f'{chrfilter_dir}/chunk_{start}_{stop}.txt'
    filteredprefix = chrfilter_dir+f'/chunk_{start}_{stop}'
    filteredprefixwithchr = chrfilter_dir+f'/chunk_{start}_{stop}.withchr'
    snpsfile = chrfilter_dir+f'/chunk_{start}_{stop}.snps'
    bcf = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.bcf'
    index = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.bcf.csi'
    chrmissing_dir = f'{qcdir}/{chrom}_missing'
    chunkmissingprefix = f'{chrmissing_dir}/chunk_{start}_{stop}'

    filteredchunk = open(chrfilterchunk).read().strip()
    print(filteredchunk)

    if filteredchunk == 'Keep':
        print(f'opening keep')

        bim = pd.read_csv(chunkbim, sep='\s+', header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2'])
        bim = bim.loc[(bim['a1'].isin(alleles+['0'])) & (bim['a2'].isin(alleles+['0']))]
        pos = np.array(bim['pos'])

        dis = pos[1:] - pos[:-1]
        if (dis < 0).sum() > 0:
            print(f'Positions are not inc in bim')
            return

        # filter A T C G
        # filter duplicates (and pos - 1 == pos)
        snp = np.array(bim['snp'])
        duplicates = np.append(snp[:-1][dis<=1], snp[1:][dis<=1])
        bim = bim.loc[~bim['snp'].isin(duplicates)]

        # filter missing rate / call rate
        lmissfile = chunkmissingprefix+'.lmiss'
        highmissing = pd.read_csv(lmissfile, sep='\s+')
        highmissing = highmissing.loc[highmissing['F_MISS'] > (1-minsnpCallRate)]
        bim = bim.loc[~bim['snp'].isin(highmissing['SNP'])]

        # extract snps
        print(snpsfile)
        bim['snp'].to_csv(snpsfile, sep='\t', index=False, header=False)
        bim['ref'] = bim['snp'].str.split(pat='_', expand=True).iloc[:,2]
        bim[['snp', 'ref']].to_csv(snpsfile+'.ref', sep='\t', index=False, header=False)

        plink_cmd = f'{PLINK} --bfile {chunkprefix} --extract {snpsfile} --maf {min_maf} --memory {mem_per_worker} --a2-allele {snpsfile}.ref 2 1 --recode vcf-iid bgz --out {filteredprefix}'
        print(plink_cmd)
        subprocess.call(plink_cmd, shell=True)
        #parse vcf to add "chr"
        with gzip.open(f'{filteredprefix}.vcf.gz', 'rt') as f:
            with gzip.open(f'{filteredprefixwithchr}.vcf.gz', 'wt') as g:
                for rline in f:
                    line = rline.strip()
                    if line.startswith('#'):
                        chr_num = chrom.replace('chr','')
                        g.write(f'{line}\n'.replace(f'ID={chr_num}', f'ID={chrom}'))
                    else:
                        g.write(f'chr{line}\n')

        subprocess.call(f'bcftools view -Ob -o {bcf} {filteredprefixwithchr}.vcf.gz', shell=True)
        subprocess.call(f'bcftools index {bcf}', shell=True)
    else:
        subprocess.call(f'touch {bcf} {index}', shell=True)


def phase_chunks(window):
    start, stop = window
    chunkprefix = chr_dir+f'/chunk_{start}_{stop}'
    chunkbim = chr_dir+f'/chunk_{start}_{stop}.bim'
    chrfilter_dir = f'{qcdir}/{chrom}_filteredChunk'
    chrfilterchunk = f'{chrfilter_dir}/chunk_{start}_{stop}.txt'
    filteredprefix = chrfilter_dir+f'/chunk_{start}_{stop}'
    bcf = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.bcf'
    index = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.bcf.csi'

    phasedprefix = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.phased'
    phasedvcf = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.phased.vcf.gz'

    filteredchunk = open(chrfilterchunk).read().strip()

    if filteredchunk == 'Keep':
        subprocess.call(f'{EAGLE} --vcfRef {ref} \
                                  --vcfTarget {bcf} \
                                  --geneticMapFile {GMAP} \
                                  --chrom {chrom} \
                                  --bpStart {start} \
                                  --bpEnd {stop} \
                                  --outPrefix {phasedprefix} \
                                  --numThreads {multicore} \
                                  --allowRefAltSwap \
                                  --vcfOutFormat z', shell=True, stdout=subprocess.DEVNULL)
    else:
        subprocess.call(f'touch {phasedvcf}', shell=True)


def impute_chunks(window):
    #done after phasing
    start, stop = window
    chunkprefix = chr_dir+f'/chunk_{start}_{stop}'
    chunkbim = chr_dir+f'/chunk_{start}_{stop}.bim'
    chrfilter_dir = f'{qcdir}/{chrom}_filteredChunk'
    chrfilterchunk = f'{chrfilter_dir}/chunk_{start}_{stop}.txt'
    filteredprefix = chrfilter_dir+f'/chunk_{start}_{stop}'
    bcf = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.bcf'
    index = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.bcf.csi'

    phasedvcf = chrfilter_dir+f'/chunk_{start}_{stop}.withchr.phased.vcf.gz'

    filteredchunk = open(chrfilterchunk).read().strip()
    buffersize = 2.5e6
    impute_dir = outdir+f'/imputation/{chrom}'
    Path(impute_dir).mkdir(parents=True, exist_ok=True)
    imputeddosefile = impute_dir+f'/chunk_{start}_{stop}.imputed.dose.vcf.gz'
    imputedprefix = impute_dir+f'/chunk_{start}_{stop}.imputed'



    if filteredchunk == 'Keep' and stop - buffersize > start + buffersize:
        # check whether variants exist in imputation region:
        # chunk 1:120000001-145000000 pass QC filter, but don't have variants in
        # imputation region 1:122000001-143000000
        empty_imputation_region = True
        for line in gzip.open(phasedvcf, 'rt'):
            if not re.search('^#', line):
                line = line.split('\t', 2)
                pos = int(line[1])
                if pos > start + buffersize and pos < stop - buffersize:
                    empty_imputation_region = False
                    break
        if empty_imputation_region:
            subprocess.call(f'touch {imputeddosefile}', shell=True)
        else:
            subprocess.call(f'{MM4} --refHaps {refm3} \
                                    --haps {phasedvcf} --chr {chrom} \
                                    --start {start} \
                                    --end {stop} \
                                    --window {windowsize} \
                                    --cpus {multicore} \
                                    --format GT,DS,GP \
                                    --allTypedSites \
                                    --meta \
                                    --noPhoneHome \
                                    --minRatio 0.00001 \
                                    --prefix {imputedprefix}', shell=True, stdout=subprocess.DEVNULL)
    else:
        subprocess.call(f'touch {imputeddosefile}', shell=True)



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

    if os.path.isfile(pos):
        pass
    else:
        #write to temp file and move if completed/successful
        subprocess.run(f"cat {tar}.bim | cut -f4 > {temp_pos}", shell=True)
        subprocess.run(f"mv {temp_pos} {pos}", shell=True)

    with open(pos) as f:
        positions = [int(line) for line in f.read().splitlines()]



    #qc target data
    qcdir = outdir+'/qc'
    Path(qcdir).mkdir(parents=True, exist_ok=True)



    #extract target chunks
    chr_dir = f'{qcdir}/{chrom}'
    if os.path.isdir(chr_dir):
        shutil.rmtree(chr_dir)
    Path(chr_dir).mkdir(parents=True, exist_ok=True)

    bim = pd.read_csv(tar+'.bim', sep='\s+', header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2'])
    windows = position_windows(pos=np.array(bim['pos']), size=int(chunksize), start=1, step=int(chunksize-overlap))

    refbim_pd = pd.read_csv(refbim, sep='\s+', header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2'])



    #with 1 worker each
    print(f'creating chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(create_chunks, windows)
        executor.shutdown()
    chunk_glob = glob.glob(chr_dir+f'/chunk_*_*.bim')
    print(f'expected: {len(windows)}, found: {len(chunk_glob)}')



    print(f'qc chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(qc_chunks, windows)
        executor.shutdown()
    chrfilter_dir = f'{qcdir}/{chrom}_filteredChunk'
    qc_glob = glob.glob(f'{chrfilter_dir}/chunk_*_*.txt')
    print(f'expected: {len(windows)}, found: {len(qc_glob)}')



    print(f'bcf chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(filter_chunks, windows)
        executor.shutdown()
    bcf_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.bcf')
    print(f'expected: {len(windows)}, found: {len(bcf_glob)}')



    #with multicore workers each
    print(f'phase chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
        executor.map(phase_chunks, windows)
        executor.shutdown()
    phase_glob = glob.glob(chrfilter_dir+f'/chunk_*_*.phased.vcf.gz')
    print(f'expected: {len(windows)}, found: {len(phase_glob)}')



    #with multicore workers each
    imputdir = outdir+'/imputation'
    Path(imputdir).mkdir(parents=True, exist_ok=True)
    print(f'impute chunks')
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_multithread) as executor:
        executor.map(impute_chunks, windows)
        executor.shutdown()
    impute_glob = glob.glob(f'{outdir}/imputation/{chrom}/chunk_*_*.imputed.dose.vcf.gz')
    print(f'expected: {len(windows)}, found: {len(impute_glob)}')

    assert len(windows) == len(impute_glob), 'num imputed chunks (%s) != num windows (%s) ' % (len(windows), len(impute_glob))



    #stitch imputed chunks together
    topmed_imputation_aggregate_chunks = [f'{outdir}/imputation/{chrom}/chunk_{start}_{stop}.imputed.dose.vcf.gz' for start,stop in windows]

    outchunks = f'{imputdir}/{chrom}.imputed.dose.vcf.chunks'

    dtype = [('start', int), ('file', 'U100')]
    chunk_vcfs = []
    for f in topmed_imputation_aggregate_chunks:
        print(f'parsing input: {f} ({os.stat(f).st_size})')
        if os.stat(f).st_size != 0:
            start = int(f.split('/')[-1].split('_')[1])
            chunk_vcfs.append((start, f))
    chunk_vcfs = np.array(chunk_vcfs, dtype=dtype)
    chunk_vcfs = np.sort(chunk_vcfs, order='start')
    chunk_vcfs = list(chunk_vcfs['file'])
    for chunk in chunk_vcfs:
        print(chunk)
    print(f'writing to {outchunks}')
    with open(outchunks, 'w') as f:
        f.write('\n'.join(chunk_vcfs))
    print(f'writing done')

    imputchromvcf = f'{imputdir}/{chrom}.imputed.dose.vcf.gz'
    subprocess.call(f'bcftools concat --threads {workers} -Oz -o {imputchromvcf} -f {outchunks}', shell=True)


