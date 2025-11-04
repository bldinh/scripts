library(argparser)

# clumpdir <- 'vanillaasegrm4anc'
# pheno <- 'bmi'

p <- arg_parser('Calc effective sample size for difference in summary stats')
p <- add_argument(p, "--dir", help="Directory containing the folders for each phenotype")
p <- add_argument(p, "--pheno", help="pheno/folder name")

argv <- parse_args(p)

clumpdir <- argv$dir
pheno <- argv$pheno

cat(paste0('clumpdir:',clumpdir,'\n'))
cat(paste0('pheno:',pheno,'\n'))

#needs at least 35gb for the columnanalysis

cutoffs <- c(1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8)

#eGRM column-wide enrichment calculation
#cutoffs

#TODO
#outfile to write to?


#grab corresponding pval
parse_pheno_by_fp <- function() {

    #read for each chrom the overall..
    cat('parsing chr1..\n')
    pvals <- read.delim(paste0(clumpdir,'/',pheno,'/',pheno,'.all.mac10.intersect.subset.bestpval.chr',1,'.plinkassocfmt.txt'), header=T, sep='\t')
    #head(pvals)
    for (i in 2:22){
       cat(paste0('chr',i,'\n'))
       next_df <- read.delim(paste0(clumpdir,'/',pheno,'/',pheno,'.all.mac10.intersect.subset.bestpval.chr',i,'.plinkassocfmt.txt'), header=T, sep='\t')
       pvals <- rbind(pvals,next_df)
    }
    cat('done!\n')

    #set datatypes
    #P	SNP	CHR	BP
    pvals$P <- as.numeric(pvals$P)
    #pvals$CHR <- as.integer(pvals$CHR)
    pvals$BP <- as.integer(pvals$BP)
    pvals$new.pval <- as.numeric(pvals$newpval)
    pvals$base.pval <- as.numeric(pvals$basepval)
    pvals <- pvals[complete.cases(pvals),]

                #orig
                ##P SNP CHR BP
                #new_pvals <- read.delim(newfp, header=T, sep=' ')
                #colnames(new_pvals) <- c('Score.pval','onpid','chr','pos','new.pval','base.pval')
                #base_pvals <- read.delim(basefp, header=T, sep=' ')
                #colnames(base_pvals) <- c('Score.pval','snpid','chr','pos','new.pval','base.pval')

                ##base fp has datatype issues when parsing
                #base_pvals$Score.pval <- as.numeric(base_pvals$Score.pval)
                #base_pvals$chr <- as.integer(base_pvals$chr)
                #base_pvals$pos <- as.integer(base_pvals$pos)
                #base_pvals <- base_pvals[complete.cases(base_pvals),]

                ##repeat for new fp
                #new_pvals$Score.pval <- as.numeric(new_pvals$Score.pval)
                #new_pvals$chr <- as.integer(new_pvals$chr)
                #new_pvals$pos <- as.integer(new_pvals$pos)
                #new_pvals <- new_pvals[complete.cases(new_pvals),]

    #parsing clumped output
    temp = list.files(path = paste0(clumpdir,"/",pheno,"/"), pattern="\\.clumped$", full.names=T)
    myfiles = lapply(temp, read.delim, sep = "")

    new_tophit_pvals <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(new_tophit_pvals) <- c('Score.pval', 'chrpos')
    base_tophit_pvals <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(base_tophit_pvals) <- c('Score.pval', 'chrpos')

    for (i in myfiles){
        #if multiple clumps
        for (rowidx in 1:nrow(i)){
            #CHR  F   SNP             BP          P       TOTAL  NSIG   S05    S01   S001  S0001  SP2
            #13   1   rs370336408     55300080    1e-11   2      0      0      0     0     2      rs140377079(1),rs142675376(1)
            chr <- i[rowidx,c("CHR")]
            pos <- i[rowidx,c("BP")]
            cur_tophits <- c(i[rowidx,'SNP'])

            SP2 <- i[rowidx,'SP2']
            if (SP2 != 'NONE'){
              poten_tophits <- stringr::str_split(SP2, ',', n = Inf, simplify = T)
              poten_tophits <- gsub("\\(.\\)", "", poten_tophits)
              for (sp2idx in 1:length(poten_tophits)){
                cur_tophits <- c(poten_tophits[sp2idx],cur_tophits)
              }
            }

            #check new
            loci_pvals <- pvals[pvals$SNP %in% cur_tophits,]
            #print(loci_pvals)
            #loci_new_pvals <- loci_new_pvals[complete.cases(loci_new_pvals),]
            #if (nrow(loci_new_pvals > 0)) {
            bestnewidx <- which.min(loci_pvals[,'new.pval'])
            bestnewpval <- loci_pvals[bestnewidx,'new.pval']
            new_tophit_pvals[nrow(new_tophit_pvals) + 1,] = c(bestnewpval, paste0(pheno, ",", chr, ",", pos))
            #}

            #check base
            #loci_base_pvals <- pvals[pvals$SNP %in% cur_tophits,'base.pval']
            #loci_base_pvals <- loci_base_pvals[complete.cases(loci_base_pvals),]
            #if (nrow(loci_base_pvals > 0)) {
            bestbaseidx <- which.min(loci_pvals[,'base.pval'])
            bestbasepval <- loci_pvals[bestbaseidx,'base.pval']
            base_tophit_pvals[nrow(base_tophit_pvals) + 1,] = c(bestbasepval, paste0(pheno, ",", chr, ",", pos))
            #}
        }
    }

    #end clumping section

    #update tophit dfs to have uniform col names
    colnames(new_tophit_pvals) <- c('new.pval','snpid')
    colnames(base_tophit_pvals) <- c('base.pval','snpid')
    tophit_pvals <- merge(new_tophit_pvals, base_tophit_pvals, by="snpid")



    #parsing top_hits
    pheno_pvals <- c()
    pheno_counts <- c()
    for (cutoff in cutoffs){

        top_hits <- union(tophit_pvals[as.numeric(tophit_pvals$new.pval) < cutoff,'snpid'], tophit_pvals[as.numeric(tophit_pvals$base.pval) < cutoff,'snpid'])
        num_snps <- length(top_hits)

        #subset
        new_tophit_pvals <- tophit_pvals[tophit_pvals$snpid %in% top_hits,'new.pval']
        base_tophit_pvals <- tophit_pvals[tophit_pvals$snpid %in% top_hits,'base.pval']

        new_z <- qnorm(as.numeric(as.character(unlist(new_tophit_pvals)))/2)
        base_z <- qnorm(as.numeric(as.character(unlist(base_tophit_pvals)))/2)

        new_chisq <- mean(new_z^2)
        base_chisq <- mean(base_z^2)

        samp_size <- (new_chisq-1)/(base_chisq-1)

        pheno_pvals <- c(pheno_pvals, samp_size)
        pheno_counts <- c(pheno_counts, length(top_hits))
    }

    cat(pheno_pvals)
    cat('\n')
    cat(pheno_counts)
    cat('\n')
    cat('\n')

    return(tophit_pvals)
}


#calc_sampsize_change <- function(basedir, newdir, phenos = c("bmi", "crp", "glucose", "hdl", "height", "hip", "homair", "insulin", "ldl", "totalcholesterol", "triglycerides", "waist", "waisthip")) {
calc_sampsize_change <- function() {

  cur_pheno_tophit_pvals <- parse_pheno_by_fp()

  #make aggregate list
  all_tophit_pvals <- cur_pheno_tophit_pvals

  #if (length(phenos) > 1) {
  #  for (x in phenos[2:length(phenos)]){
  #    cat(x)
  #    cat('\n')
  #    newfp <- paste0(newdir,"/",x,"/allchr.fullvariant.pval.snpid.chr.pos.txt")
  #    basefp <- paste0(basedir,"/",x,"/allchr.fullvariant.pval.snpid.chr.pos.txt")
  #    cur_pheno_tophit_pvals <- parse_pheno_by_fp(basefp, basedir, newfp, newdir, x)
  #    all_tophit_pvals <- rbind(all_tophit_pvals, cur_pheno_tophit_pvals)
  #  }
  #}

  cat('overall:\n')
  pheno_pvals <- c()
  pheno_counts <- c()
  all_tophit_pvals$new.pval <- as.numeric(all_tophit_pvals$new.pval)
  all_tophit_pvals$base.pval <- as.numeric(all_tophit_pvals$base.pval)
  #all_tophit_pvals_complete <- all_tophit_pvals[complete.cases(all_tophit_pvals),]
  all_tophit_pvals_complete <- all_tophit_pvals

  for (cutoff in cutoffs){
   #print(cutoff)
   top_hits <- union(all_tophit_pvals_complete[as.numeric(all_tophit_pvals_complete$new.pval) < cutoff, 'snpid'], 
                     all_tophit_pvals_complete[as.numeric(all_tophit_pvals_complete$base.pval) < cutoff, 'snpid'])
   num_snps <- length(top_hits)

   #subset
   new_tophit_pvals <- all_tophit_pvals_complete[all_tophit_pvals_complete$snpid %in% top_hits,'new.pval']
   base_tophit_pvals <- all_tophit_pvals_complete[all_tophit_pvals_complete$snpid %in% top_hits,'base.pval']

   new_z <- qnorm(as.numeric(as.character(unlist(new_tophit_pvals)))/2)
   base_z <- qnorm(as.numeric(as.character(unlist(base_tophit_pvals)))/2)

   new_chisq <- mean(new_z^2)
   base_chisq <- mean(base_z^2)

   samp_size <- (new_chisq-1)/(base_chisq-1)
   pheno_pvals <- c(pheno_pvals, samp_size)
   pheno_counts <- c(pheno_counts, num_snps)
  }

  #print(all_tophit_pvals_complete[all_tophit_pvals_complete$snpid %in% top_hits,])

  cat(pheno_pvals)
  cat('\n')
  cat(pheno_counts)
  cat('\n')
  cat('\n')
  ms <- gc()
  cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
  
  saveRDS(all_tophit_pvals_complete, file = paste0(clumpdir,"/",pheno,"/all_tophit_pvals_complete.rds"))
  outfile <- paste0(clumpdir,'/',pheno,'/','samplesize.txt')
  write.table(as.list(cutoffs), outfile, sep='\t', row.names = F, col.names = F)
  formatted_pheno_pvals <- mapply(round, pheno_pvals, digits=5)
  write.table(as.list(formatted_pheno_pvals), outfile, sep='\t', row.names = F, col.names = F, append = T)
  write.table(as.list(pheno_counts), outfile, sep='\t', row.names = F, col.names = F, append = T)
}



#calc_sampsize_change('vanilla','asegrm4anc', phenos=c('bmi'))
calc_sampsize_change()

