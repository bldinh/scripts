#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --time=7-00:00:00
#SBATCH --array=1-22
#SBATCH --cpus-per-task=10
#SBATCH --mem=20gb

# Takes in phased vcf to be used as reference and:
#   Runs tabix on vcf
#   Creates m3vcf using Minimac3
#   Create bim (only) using PLINK

threads=10
mem=20000

PLINK=path/to/plink
MINIMAC3=path/to/mm3


vcf=path/to/phased/and/gzipped/ref/vcf
outdir=path/to/dir/to/write/to
outprefix=somePrefix

#
# changes below this line are not required
#

mkdir -p $outdir
mkdir -p $outdir/ref/vcf
mkdir -p $outdir/ref/bim
mkdir -p $outdir/ref/m3vcfs

refprefix=$outdir/ref/vcf/chr$SLURM_ARRAY_TASK_ID.$outprefix.phased
refphasedvcfgz=$refprefix.vcf.gz

ln -s $vcf $refphasedvcfgz
tabix -p vcf $refphasedvcfgz

m3vcfprefix=$outdir/ref/m3vcfs/chr$SLURM_ARRAY_TASK_ID.$outprefix.phased
$MINIMAC3 --refHaps $refphasedvcfgz --processReference --prefix $m3vcfprefix --myChromosome chr$SLURM_ARRAY_TASK_ID

$PLINK --double-id \
       --keep-allele-order \
       --make-just-bim \
       --memory $mem \
       --out $outdir/ref/bim/chr$SLURM_ARRAY_TASK_ID.$outprefix.phased \
       --vcf $refphasedvcfgz \
       --vcf-half-call missing


