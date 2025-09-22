library("argparser")

options(scipen=10000)
bg="white"

parser <- arg_parser("")

parser <- add_argument(parser,"--ref", help="the reference chr pos file")
parser <- add_argument(parser,"--tar", help="the target chr pos file")
parser <- add_argument(parser,"--out", help="output filename")
parser <- add_argument(parser,"--prefix", help="plot title")

argv <- parse_args(parser)

ref_fp <- argv$ref
tar_fp <- argv$tar
prefix <- argv$prefix

ref_pos <- read.csv(ref_fp, header=F, sep='\t')
tar_pos <- read.csv(tar_fp, header=F, sep='\t')

last_pos <- max(ref_pos$V1,tar_pos$V1)
max_mb <- (last_pos %/% 1000000) + 1

pdf(argv$out)

p1 <- hist(ref_pos$V1, breaks=max_mb, col=rgb(1,0,0,.75), main=paste(prefix, '(', length(ref_pos$V1)-length(tar_pos$V1), 'variants lost)'))
p2 <- hist(tar_pos$V1, breaks=max_mb, col=rgb(0,0,1,.25), add=T)

ref_title <- paste('ref (', length(ref_pos$V1), 'variants )')
tar_title <- paste('tar (', length(tar_pos$V1), 'variants )')

legend("bottomleft", inset=.02, legend=c(ref_title, tar_title), fill=c('red','blue'), cex=0.8)

dev.off()




