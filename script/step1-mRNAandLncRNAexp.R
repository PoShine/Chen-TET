rm(list = ls())
library(dplyr)
##读取RNA-seq表达矩阵，htseq_log2(fpkm+1)
setwd("F:/chenkg/THYM/data")
htseq_fpkm = read.table("TCGA-THYM.htseq_fpkm.tsv",header = T,row.names = 1,
                        sep = "\t",stringsAsFactors = F,check.names = F)
##读取RNA-seq表达矩阵，htseq_log2(count+1)
htseq_counts = read.table("TCGA-THYM.htseq_counts.tsv",header = T,row.names = 1,
                        sep = "\t",stringsAsFactors = F,check.names = F)
##读取gencode.v22.annotation.gene注释文件
gencode.v22 = read.table("1gencode.v22.annotation.gene.probeMap",header = T,
                         sep = "\t",stringsAsFactors = F)

##extract mRNA expression
if(T){
  gencode.v22.mRNA = filter(gencode.v22,genetype=="protein_coding")[,1:2]
  #fpkm
  mRNA_exp_fpkm = htseq_fpkm[gencode.v22.mRNA$id,]
  mRNA_exp_fpkm$id = rownames(mRNA_exp_fpkm)
  mRNA_exp_fpkm = left_join(mRNA_exp_fpkm,gencode.v22.mRNA,by="id")
  mRNA_expr_fpkm = as.data.frame(t(sapply(split(mRNA_exp_fpkm,mRNA_exp_fpkm$gene),function(x) colMeans(x[,1:(ncol(x)-2)]))),stringsAsFactors = F)
  #log2(count+1) to count
  mRNA_exp_counts = htseq_counts[gencode.v22.mRNA$id,]
  mRNA_exp_counts1 = 2^mRNA_exp_counts-1
  mRNA_exp_counts1$id = rownames(mRNA_exp_counts1)
  mRNA_exp_counts1 = left_join(mRNA_exp_counts1,gencode.v22.mRNA,by="id")
  mRNA_exp_counts = as.data.frame(t(sapply(split(mRNA_exp_counts1,mRNA_exp_counts1$gene),function(x) colMeans(x[,1:(ncol(x)-2)]))),stringsAsFactors = F)
}


##extract lncRNA expression
if(T){
  gencode.v22.lncRNA = filter(gencode.v22,
                              genetype=="lincRNA" |
                                genetype=="antisense" |
                                genetype=="sense_intronic" |
                                genetype=="3prime_overlapping_ncrna" |
                                genetype=="sense_overlapping" |
                                genetype=="non_coding" |
                                genetype=="macro_lncRNA")[,1:2]
  #fpkm
  lncRNA_exp_fpkm = htseq_fpkm[gencode.v22.lncRNA$id,]
  lncRNA_exp_fpkm$id = rownames(lncRNA_exp_fpkm)
  lncRNA_exp_fpkm = left_join(lncRNA_exp_fpkm,gencode.v22.lncRNA,by="id")
  lncRNA_exp_fpkm = as.data.frame(t(sapply(split(lncRNA_exp_fpkm,lncRNA_exp_fpkm$gene),function(x) colMeans(x[,1:(ncol(x)-2)]))),stringsAsFactors = F)
  #log2(count+1) to count
  lncRNA_exp_counts = htseq_counts[gencode.v22.lncRNA$id,]
  lncRNA_exp_counts1 = 2^lncRNA_exp_counts-1
  lncRNA_exp_counts1$id = rownames(lncRNA_exp_counts1)
  lncRNA_exp_counts1 = left_join(lncRNA_exp_counts1,gencode.v22.lncRNA,by="id")
  lncRNA_exp_counts = as.data.frame(t(sapply(split(lncRNA_exp_counts1,lncRNA_exp_counts1$gene),function(x) colMeans(x[,1:(ncol(x)-2)]))),stringsAsFactors = F)
}

##保存mRNA_expr_fpkm、mRNA_exp_counts、lncRNA_exp_fpkm、lncRNA_exp_counts
save(mRNA_expr_fpkm,mRNA_exp_counts,lncRNA_exp_fpkm,lncRNA_exp_counts,file = "mRNAandLncRNAexp.RData")

##读取miRNA表达数据log2(RPM+1)
miRNA_exp_RPM = read.table("TCGA-THYM.mirna.tsv",header = T,row.names = 1,
                           sep = "\t",stringsAsFactors = F,check.names = F)
##读取miRNA表达数据count
miRNA_exp_count = read.csv("TCGA-THYM-miRNA.csv",header = T,row.names = 1,stringsAsFactors = F,check.names = F)
colnames(miRNA_exp_count) = sub("read_count_","",colnames(miRNA_exp_count))
colnames(miRNA_exp_count) = sub("-[123][12]R-.*","",colnames(miRNA_exp_count))
miRNA_exp_count = miRNA_exp_count[,colnames(miRNA_exp_RPM)]

save(mRNA_expr_fpkm,mRNA_exp_counts,lncRNA_exp_fpkm,lncRNA_exp_counts,miRNA_exp_count,miRNA_exp_RPM,file = "mRNAandLncRNAandmiRNAexp.RData")

