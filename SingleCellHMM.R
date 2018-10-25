# R --vanilla --slave --args $(pwd) PREFIX.bed < SingleCellHMM.R 

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_f = args[2] #PREFIX.bed.gz #05-007B1_gene_exon_tagged.REF_chr22_split.bed.gz



library(groHMM)
library(rtracklayer)

GRangeTobed <- function(gr, f_name){
  df <- data.frame(seqnames=seqnames(gr),
    starts=start(gr)-1,
    ends=end(gr),
    names=c(rep(".", length(gr))),
    scores=c(rep(".", length(gr))),
    strands=strand(gr))

  write.table(df, file=f_name, quote=F, sep="\t", row.names=F, col.names=F)
}


S_split =  import(input_f)
hmmResult_AllEM_split <- detectTranscripts_AllEM(S_split)
GRangeTobed(hmmResult_AllEM_split$transcripts, paste(strsplit(input_f, '.bed', T)[[1]][1], "_HMM.bed", sep = ""))

