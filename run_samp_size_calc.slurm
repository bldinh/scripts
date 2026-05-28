#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=0-12
#SBATCH --time=0-10:00:00
#SBATCH --output=logs/samplesize/pheno%a.out
#SBATCH --error=logs/samplesize/pheno%a.out
#SBATCH --job-name=sampsize
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --partition=

module load r/4.4.2

#        0     1     2         3     4        5     6        7         8     9                  10              11      12         13             14
myArray=("bmi" "crp" "glucose" "hdl" "height" "hip" "homair" "insulin" "ldl" "totalcholesterol" "triglycerides" "waist" "waisthip" "waisthipsex1" "waisthipsex2")
pheno=${myArray[$SLURM_ARRAY_TASK_ID]}

outputdir=vanillaeGRMwithPCs

#Summary stats expected to have the following fields:
    #snpid
    #Score.pval
    #Rsq
    #freq
    #MAC
    #chr
    #pos
#if column names differ, they can be provided to the first script
#default filtering is at MAF >= 0.005, MAC >= 10, Rsq >= 0.3

#if all phenotypes are in the same directory, the paths below can use arrays to run multiple phenotypes at once

new_assoc=/path/to/assoc/sum/stats/chr%s.addRsq.addAF.txt
base_assoc=/path/to/baseline/assoc/sum/stats/chr%s.addRsq.addAF.txt


python3 _1_calc_col_enrich_from_sumstat.py --pheno $pheno \
                                           --path1 $new_assoc \
                                           --path2 $base_assoc \
                                           --outdir $outputdir

Rscript _2_calc_eff_samp_size_ldpruned.R --dir $outputdir --pheno $pheno


