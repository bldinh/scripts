library("argparser")

options(scipen=10000)
bg="white"

parser <- arg_parser("")

parser <- add_argument(parser,"--ref", help="the reference chr pos file")
parser <- add_argument(parser,"--tar", help="the target chr pos file")
parser <- add_argument(parser,"--out", help="output filename")
parser <- add_argument(parser,"--prefix", help="plot title")

argv <- parse_args(parser)

ref_pos <- read.csv(argv$ref, header=F, sep='\t')
tar_pos <- read.csv(argv$tar, header=F, sep='\t')

last_pos <- max(ref_pos$V1,tar_pos$V1)
max_mb <- (last_pos %/% 1000000) + 1
#hist_breaks <- c(1,c(1:max_mb)*1000000)

pdf(argv$out)

p1 <- hist(ref_pos$V1, breaks=max_mb, col=rgb(1,0,0,.75), main=argv$prefix)
p2 <- hist(tar_pos$V1, breaks=max_mb, col=rgb(1,0.5,0,.75), add=T)

legend("bottomleft", inset=.02, legend=c('ref', 'tar'), fill=c('red','orange'), cex=0.8)

dev.off()




