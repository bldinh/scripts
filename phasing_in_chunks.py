import argparse
from concurrent.futures import ThreadPoolExecutor
import glob
import numpy as np
import os
from pathlib import Path
import shutil
import subprocess


# based on implementation of the TOPMed pipeline by Minhui Chen


#
# Program paths
#


EAGLE = '/scratch1/bldinh/programs/eagle'


#
# Parse arguments
#


parser = argparse.ArgumentParser()
parser.add_argument('--tar', required=True, help='target data in vcf format, chrom format should match ref. chrom format (if phasing vs. ref.)')
parser.add_argument('--chrom', required=True, help='should be formatted the way "chrom" is in target vcf, e.g. "chr22" or "22"')
parser.add_argument('--outdir', required=True, help='without trailing "/", will be made if it does not exist')
parser.add_argument('--workers', required=True, help='number of chunks to work on at same time')
parser.add_argument('--ref', required=False, help='phased reference vcf', default='')
parser.add_argument('--chunksize', required=False, help='chunk size in bp', default=25e6)
parser.add_argument('--overlap', required=False, help='overlap size in bp', default=5e6)
parser.add_argument('--multicore', required=False, help='number of workers to use when multithreading', default=5)
parser.add_argument('--gmap', required=False, help='genetic map to run Eagle with', default='/scratch1/bldinh/programs/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz')

args = parser.parse_args()
tar = args.tar
OUTDIR = args.outdir
chrom = args.chrom
GMAP = args.gmap
ref = args.ref
workers = int(args.workers)
chunksize = int(args.chunksize)
overlap = int(args.overlap)
multicore = int(args.multicore)
max_multithread = max(1,int(workers/multicore))


#
# Functions
#


def position_windows(pos, size, start=None, stop=None, step=None):
    """
    Take in list of positions and break into intervals
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


def create_chunks_and_idx(window):
    """
    Take in vcf and:
        Use bcftools to extract region between start and end and output as bcf
        Index bcf file
    """
    start, stop = window
    bcf = CHRDIR+f'/chunk_{start}_{stop}.bcf'
    subprocess.call(f'bcftools view -t {chrom}:{start}-{stop} -Ob -o {bcf} {tar}', shell=True, stdout=subprocess.DEVNULL)
    subprocess.call(f'bcftools index {bcf}', shell=True)


def phase_chunks(window):
    """
    Phase bcf chunk with Eagle and output as compressed vcf
    """
    start, stop = window
    prefix = CHRDIR+f'/chunk_{start}_{stop}'
    bcf = f'{prefix}.bcf'
    phasedprefix = f'{prefix}.phased'

    if len(ref) > 0:
        ref_tar_arg = f' --vcfRef {ref} --vcfTarget {bcf}'
    else:
        ref_tar_arg = f'--vcf {bcf}'
    subprocess.call(f'{EAGLE} {ref_tar_arg} \
                              --geneticMapFile {GMAP} \
                              --chrom {chrom} \
                              --bpStart {start} \
                              --bpEnd {stop} \
                              --outPrefix {phasedprefix} \
                              --numThreads {multicore} \
                              --allowRefAltSwap \
                              --vcfOutFormat z', shell=True, stdout=subprocess.DEVNULL)


def resize_chunks(window_tuple):
    """
    Prepare phased vcf for stitching by resizing/trimming regions
    """
    window, new_window = window_tuple
    start, stop = window
    new_start, new_stop = new_window

    phasedprefix = CHRDIR+f'/chunk_{start}_{stop}.phased'
    vcf = f'{phasedprefix}.vcf.gz'
    new_vcf = f'{phasedprefix}.resizeBin.vcf.gz'

    subprocess.call(f'bcftools view -t {chrom}:{new_start}-{new_stop} -Oz -o {new_vcf} {vcf}', shell=True)


#
# Start
#


if __name__ == '__main__':


    # create outfolder
    Path(OUTDIR).mkdir(parents=True, exist_ok=True)


    #
    # Create pos file
    #


    # create "pos" folder
    POSDIR = OUTDIR+'/targetpositions'
    Path(POSDIR).mkdir(parents=True, exist_ok=True)


    # create pos file
    pos = POSDIR+f'/{chrom}.pos.txt'
    temp_pos = pos+'.temp'

    #write to temp file and move if completed/successful
    with open(temp_pos, 'w') as f:
        if tar.endswith('.gz'):
            subprocess.run(f"zcat {tar} | grep ^[^#] | cut -f2", shell=True, stdout=f)
        else:
            subprocess.run(f"cat {tar} | grep ^[^#] | cut -f2", shell=True, stdout=f)
    subprocess.run(f"mv {temp_pos} {pos}", shell=True)

    with open(pos) as f:
        positions = [int(line) for line in f.read().splitlines()]
    windows = position_windows(pos=np.array(positions), size=int(chunksize), start=1, step=int(chunksize-overlap))


    #
    # Create chunks
    #


    # create chunk dir
    CHUNKDIR = OUTDIR+'/unphased_chunks'
    Path(CHUNKDIR).mkdir(parents=True, exist_ok=True)


    # create subfolder in chunk dir for current chromosome
    CHRDIR = f'{CHUNKDIR}/{chrom}'
    if os.path.isdir(CHRDIR):
        shutil.rmtree(CHRDIR)
    Path(CHRDIR).mkdir(parents=True, exist_ok=True)


    # extract chunks in parallel (tied to num. of workers provided)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(create_chunks_and_idx, windows)
        executor.shutdown()
    chunk_glob = glob.glob(CHRDIR+f'/chunk_*_*.bcf')
    print(f'expected: {len(windows)}, found: {len(chunk_glob)}')
    assert len(windows) == len(chunk_glob), 'num. of chunks (%s) does not match expected num. (%s)' % (len(chunk_glob), len(windows))


    #
    # Phase chunks
    #


    # phase chunks in parallel (more workers per chunk)
    with ThreadPoolExecutor(max_workers=max_multithread) as executor:
        executor.map(phase_chunks, windows)
        executor.shutdown()
    phase_glob = glob.glob(CHRDIR+f'/chunk_*_*.phased.vcf.gz')
    print(f'expected: {len(windows)}, found: {len(phase_glob)}')
    assert len(windows) == len(phase_glob), 'num. phased chunks (%s) does not match expected num. (%s) ' % (len(phase_glob), len(windows))


    # create output dir for phased results
    PHASEDIR = OUTDIR+'/phased'
    Path(PHASEDIR).mkdir(parents=True, exist_ok=True)


    #
    # Stitch phased chunks into one file
    #


    # collect phased chunk files

    new_windows = []
    for idx,interval in enumerate(windows):
        start,stop = interval
        if start == windows[0][0] and stop == windows[0][1]:
            stop = int(stop) -  int(overlap/2)
        elif start == windows[-1][0] and stop == windows[-1][1]:
            start = int(start) + int(overlap/2)
        else:
            start = int(start) + int(overlap/2)
            stop = int(stop) - int(overlap/2)
        new_windows.append([start,stop])


    # resize/trim phased chunks (parallel)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(resize_chunks, zip(windows, new_windows))
        executor.shutdown()
    resized_glob = glob.glob(CHRDIR+f'/chunk_*_*.phased.resizeBin.vcf.gz')
    print(f'expected: {len(new_windows)}, found: {len(resized_glob)}')
    assert len(new_windows) == len(resized_glob), 'num. of chunks (%s) does not match expected num. (%s)' % (len(resized_glob), len(new_windows))


    # sort phased chunk files
    aggregate_phased_chunks = [f'{CHRDIR}/chunk_{start}_{stop}.phased.resizeBin.vcf.gz' for start,stop in windows]
    outchunks = f'{PHASEDIR}/{chrom}.phased.vcf.chunks'
    chunk_vcfs = []
    for path in aggregate_phased_chunks:
        print(f'parsing input: {path} ({os.stat(path).st_size})')
        if os.stat(path).st_size != 0:
            start = int(path.split('/')[-1].split('_')[1])
            chunk_vcfs.append((start, path))
    dtype = [('start', int), ('file', 'U300')]
    chunk_vcfs = np.array(chunk_vcfs, dtype=dtype)
    chunk_vcfs = np.sort(chunk_vcfs, order='start')
    chunk_vcfs = list(chunk_vcfs['file'])
    for chunk in chunk_vcfs:
        print(chunk)


    # write chunks to file
    print(f'writing to {outchunks}')
    with open(outchunks, 'w') as f:
        f.write('\n'.join(chunk_vcfs))
    print(f'finished writing')


    # combine phased chunks
    phasedchromfile = f'{PHASEDIR}/{chrom}.phased.vcf.gz'
    subprocess.call(f'bcftools concat --threads {workers} -Oz -o {phasedchromfile} -f {outchunks}', shell=True)

    if os.path.isfile(phasedchromfile):
        shutil.rmtree(CHRDIR)


