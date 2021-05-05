rm(list = ls())
library(dplyr)

setwd("F:/chenkg/THYM")
mRNAs = c("CCDC80","DOCK11","FBN1","GALNT16","HAND2","MCAM","MYO10","NTRK2","PDE3B","SLC2A14","WASF3")
mRNA.interact = as.data.frame(t(combn(mRNAs,2)),stringsAsFactors = F)

ppi = read.table("protein_links_symbol.txt",header = F,sep = "\t",stringsAsFactors = F)


index = NULL
for (i in 1:dim(mRNA.interact)[1]) {
  print(i)
  a1 = filter(ppi,V1==mRNA.interact[i,1],V2==mRNA.interact[i,2])
  a2 = filter(ppi,V1==mRNA.interact[i,2],V2==mRNA.interact[i,1])
  a = rbind(a1,a2)
  if(dim(a)[1] >= 1){
    index = c(index,i)
  }
}










