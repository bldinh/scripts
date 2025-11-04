import argparse
import gzip
import os.path
import os
import subprocess


PLINK = '/project/chia657_28/programs/plink/plink1.90/plink'


parser = argparse.ArgumentParser()
parser.add_argument('--pheno')
parser.add_argument('--path1', help='Full path to first/new assoc file with %s instead of chromosome number')
parser.add_argument('--path2', help='Full path to second/baseline assoc file with %s instead of chromosome number')
parser.add_argument('--outdir', help='Output directory (a new folder for the phenotype will be made in the directory).')

#default settings for GENESIS sumstats run on TOPMed imputed dosages
parser.add_argument('--vcf', help='Path to imputed vcf to calculate LD.', default=f'/project2/chia657_28/datasets/private/MEC/Imputed_genotypes/NH_megagda_topmed_r3_imputed/merged.megaandgda.imputed.chr%s.vcf.gz')
parser.add_argument('--snpcol', default='snpid')
parser.add_argument('--pcol', default='Score.pval')
parser.add_argument('--rsqcol', default='Rsq')
parser.add_argument('--afcol', default='freq')
parser.add_argument('--maccol', default='MAC')
parser.add_argument('--chrcol', default='chr')
parser.add_argument('--poscol', default='pos')
parser.add_argument('--mac', default=10)
parser.add_argument('--rsq', default=0.3)
parser.add_argument('--maf', default=0.005)
parser.add_argument('--minpvalclump', default=0.00001)

args = parser.parse_args()
pheno = args.pheno
path1 = args.path1
path2 = args.path2
outdir = args.outdir
snpcol = args.snpcol
pcol = args.pcol
rsqcol = args.rsqcol
maccol = args.maccol
afcol = args.afcol
chrcol = args.chrcol
poscol = args.poscol
minmac = int(args.mac)
minrsq = float(args.rsq)
minmaf = float(args.maf)
minpvalclump = float(args.minpvalclump)


#previous script:
#    (1) preprocess sum stats and output a file
#        filter: MAC >= 10, Rsq >= 0.3, 0.005 < AF < 0.995
#    (2) parse output file from (1) and run PLINK clumping
#    (3) run chisq comparison at different cutoffs

#new script:
#    (1) (a) parse sum stats from both runs and filter via streaming
#        (b) compare pval genome-wide and write plink file
#    (2) run chisq comparison at different cutoffs


def parse_assoc_filter_and_return_dict(path):
    #take in assoc path with header
    #filter for AF, MAC, and Rsq
    #return dict with snpid as key

    d = {}
    d_all = {}

    with open(path) as f:
        parsed_header = False
        for rline in f:
            line = rline.split()

            if not parsed_header:
                afidx = line.index(afcol)
                chridx = line.index(chrcol)
                posidx = line.index(poscol)
                macidx = line.index(maccol)
                pidx = line.index(pcol)
                rsqidx = line.index(rsqcol)
                snpidx = line.index(snpcol)
                parsed_header = True
            else:
                cur_af = float(line[afidx])
                cur_mac = int(line[macidx])
                cur_chr = int(line[chridx])
                cur_pos = int(line[posidx])
                cur_rsq = float(line[rsqidx])
                cur_pval = float(line[pidx])
                cur_snpid = line[snpidx]


                if minmaf <= cur_af <= 1-minmaf and cur_mac >= minmac and cur_rsq >= minrsq:
                    d_all[cur_snpid] = [cur_pval, cur_snpid, cur_chr, cur_pos]

                    if cur_pval < minpvalclump:
                        d[cur_snpid] = [cur_pval, cur_snpid, cur_chr, cur_pos]
    return d, d_all


print(pheno)
failed_clump_list = []

subprocess.call(f'mkdir -p {outdir}/{pheno}'.split())

new_all_dict = {}
base_all_dict = {}

for chrom in range(1,22+1):
    #iterate over chromosomes 1-22
    print(f'chr{chrom}')


    #lookup dictionary logic to compare summary stats
    #read assoc1 into a dictionary
    new_assoc = path1 % (chrom)
    new_dict, new_dict_all = parse_assoc_filter_and_return_dict(new_assoc)
    print(f'    num lines {new_assoc} tophits: {len(new_dict.keys())}')
    print(f'    num lines {new_assoc} all: {len(new_dict_all.keys())}')

    #read assoc2 into a dictionary
    base_assoc = path2 % (chrom)
    base_dict, base_dict_all = parse_assoc_filter_and_return_dict(base_assoc)
    print(f'    num lines {base_assoc} tophits: {len(base_dict.keys())}')
    print(f'    num lines {base_assoc} all: {len(base_dict_all.keys())}')

    set_found = set(new_dict.keys()).union(set(base_dict.keys()))
    print(f'    num union tophits: {len(set_found)}')
    print(f'    num intersect all: {len(set(new_dict_all.keys()) & set(base_dict_all.keys()))}')


    #iterate over 2 dictionaries
    best_snps = [] #this is actually going to be all snps, not just the most significant ones
    for snpid in set(new_dict_all.keys()) & set(base_dict_all.keys()):
        new = new_dict_all[snpid]
        base = base_dict_all[snpid]
        new_pval = new[0]
        base_pval = base[0]
        if new_pval <= base_pval:
            best_snps.append(new + [new_pval, base_pval])
        elif base_pval < new_pval:
            best_snps.append(base + [new_pval, base_pval])

    fp = f'{outdir}/{pheno}/{pheno}.all.mac10.intersect.subset.bestpval.chr{chrom}.plinkassocfmt.txt'
    with open(fp, 'w') as g:
        g.write(f'P\tSNP\tCHR\tBP\tnewpval\tbasepval\n')
        for l in best_snps:
            l_chrom = l[2]
            l_chrom = list(l)
            l_chrom[2] = f'chr{chrom}'
            g.write('\t'.join([str(v) for v in l_chrom])+'\n')

    txtfile = f'{outdir}/{pheno}/{pheno}.1e5.mac{minmac}.maf{minmaf}.rsq{minrsq}.subset.bestpval.chr{chrom}.regionfile.txt'
    with open(txtfile, 'w') as g:
        for l in best_snps:
            if l[0] < minpvalclump:
                g.write(f'chr{chrom}\t{l[3]}\n')


    #extract from vcf
    outsubset = f'{outdir}/{pheno}/{pheno}.1e5.mac{minmac}.maf{minmaf}.rsq{minrsq}.subset.bestpval.chr{chrom}.vcf.gz'
    newoutsubset = f'{outdir}/{pheno}/{pheno}.1e5.mac{minmac}.maf{minmaf}.rsq{minrsq}.subset.bestpval.chr{chrom}.newid.vcf.gz'
    imputed_vcf = args.vcf % (chrom)

    cmd = f'bcftools view -i R2>.9 -R {txtfile} -o {outsubset} {imputed_vcf}'
    subprocess.call(cmd.split())


    #check vcf if there are any variants
    cmd = f'zcat {outsubset} | grep ^[^#] | wc -l'
    wc_output = subprocess.check_output(cmd, shell=True)
    wc_num = int(wc_output.decode().strip().split()[0])


    #run clumping
    if wc_num > 0:
        #quick hack to edit the vcf snpids
        with gzip.open(outsubset, 'rt') as f:
            with gzip.open(newoutsubset, 'wt') as g:
                prev_snpid = ''
                for rline in f:
                    if rline.startswith('#'):
                        g.write(rline)
                    else:
                        line = rline.strip().split()
                        cur_snpid = line[2]
                        if cur_snpid == prev_snpid:
                            line[2] = f'{prev_snpid}-{line[5]}'
                        g.write('\t'.join(line)+'\n')
                        prev_snpid = cur_snpid

        outprefix = f'{outdir}/{pheno}/{pheno}.1e5.mac{minmac}.maf{minmaf}.rsq{minrsq}.subset.bestpval.chr{chrom}'
        clump_file = outprefix+'.clumped'

        try:
            os.remove(clump_file)
        except OSError:
            pass

        cmd = f'{PLINK} --vcf {newoutsubset} --clump {fp} --clump-r2 0.1 --clump-p1 {minpvalclump} --out {outprefix}'
        log = subprocess.check_output(cmd.split())

        if os.path.isfile(clump_file):
            pass
        else:
            print(f'Failed to created file: {outprefix}.clumped')
            failed_clump_list.append([pheno, chrom])
    else:
        print('skipping clumping: no variants in vcf pass rsq filter')

    print()

print('\n\n\n')

if len(failed_clump_list) == 0:
    print('No failed clumping!')
else:
    print('Failed clumps:')
for failed in failed_clump_list:
    print(failed)

#calc chisq for variants and write to file



