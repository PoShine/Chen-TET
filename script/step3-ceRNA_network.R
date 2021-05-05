rm(list = ls())
library(dplyr)
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library (VennDiagram)
library(data.table)
setwd("F:/chenkg/THYM/data")

##获得miRNA靶点信息
if(F){
  ####miRNA-mRNA
  miRTarBase = read.table("miRNA_target/miRTarBase_MTI.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F,check.names = F)
  miRTarBase = filter(miRTarBase,miRTarBase$`Species (miRNA)`=="Homo sapiens")
  miRTarBase_tar = miRTarBase[,c(2,4)]
  colnames(miRTarBase_tar) = c("miRNA","mRNA")
  starbase = read.table("miRNA_target/starBaseV3_miRNA_mRNA.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F,check.names = F)
  starbase_tar = starbase[,c(2,4)]
  colnames(starbase_tar) = c("miRNA","mRNA")
  miRNA_mRNA = distinct(rbind(miRTarBase_tar,starbase_tar))
  miRNA_mRNA$miRNA = apply(miRNA_mRNA,1,
                           function(x){
                             return(str_c(str_split(x[1],"-")[[1]][1:3],collapse = "-"))
                           })
  miRNA_mRNA = distinct(miRNA_mRNA)
  
  ####miRNA-lncRNA
  if(F){
    lncBase = fread("miRNA_target/lncBaseV2_predicted_data/lncBaseV2_predicted_human_data.csv",
                    sep = "\t",quote = "")
    lncBase = lncBase[grep(">",unlist(lncBase[,1])),]
    write.table(lncBase,file = "lncbase_miRNA_lncRNA.txt",quote = F,col.names = F)
  }
  lncBase = fread("miRNA_target/lncbase_miRNA_lncRNA.txt",sep = ",",quote = "")
  lncBase = filter(lncBase, V4>=0.8)
  lncBase = filter(lncBase,grepl("ENSG",V2))
  lncBase = lncBase[,3:2]
  colnames(lncBase) = c("miRNA","lncRNA")
  lncBase$miRNA = sub("\\(.*\\)","",lncBase$miRNA)
  lncBase$lncRNA = sub("ENSG.*\\(","",lncBase$lncRNA)
  lncBase$lncRNA = sub("\\)","",lncBase$lncRNA)
  
  starbaselnc = read.table("miRNA_target/starBaseV3_miRNA_lncRNA.txt",header = T,sep = "\t",quote = "",stringsAsFactors = F,check.names = F)
  miRNA_lncRNA_star = starbaselnc[,c(2,4)]
  colnames(miRNA_lncRNA_star) = c("miRNA","lncRNA")
  miRNA_lncRNA = rbind(miRNA_lncRNA_star,lncBase)
  miRNA_lncRNA$miRNA = apply(miRNA_lncRNA,1,
                           function(x){
                             return(str_c(str_split(x[1],"-")[[1]][1:3],collapse = "-"))
                           })
  miRNA_lncRNA = distinct(miRNA_lncRNA)
  
  save(miRNA_mRNA,miRNA_lncRNA,file = "miRNA_target.RData")
}

load("miRNA_target.RData")
load("step2-DEGs.RData")

###1.mRNAs和lncRNAs是ceRNA；2.mRNAs和lncRNAs表达的相关性为正相关
diff_mRNAs
diff_lncRNAs
diff_miRNAs

##相关性计算
library(WGCNA)
diff_mRNAs_fpkm_t = mRNA_expr_fpkm[diff_mRNAs,c(colnames(mRNA_carcinoma),colnames(mRNA_typeA),colnames(mRNA_typeB),colnames(mRNA_typeAB))]
diff_lncRNAs_fpkm_t = lncRNA_exp_fpkm[diff_lncRNAs,c(colnames(lncRNA_carcinoma),colnames(lncRNA_typeA),colnames(lncRNA_typeB),colnames(lncRNA_typeAB))]
mRNA_lncRNA_cor = corAndPvalue(t(diff_mRNAs_fpkm_t),t(diff_lncRNAs_fpkm_t))

##超几何检验
mRNA_lncRNA_phyper = function(mrna,lncrna){
  mRNA_tar = filter(miRNA_mRNA,mRNA==mrna)
  lncRNA_tar = filter(miRNA_lncRNA,lncRNA==lncrna)
  N = 1881
  m = length(mRNA_tar$miRNA)
  n = N-m
  k = length(lncRNA_tar$miRNA)
  q = length(intersect(mRNA_tar$miRNA,lncRNA_tar$miRNA))
  p = phyper(q-1, m, n, k, lower.tail = F)
  return(p)
}

mRNA_ce_lncRNA = NULL
mRNAs = c()
lncRNAs = c()
phyper_p = c()
cor_coe = c()
cor_p = c()
for(m in diff_mRNAs){
  for(lnc in diff_lncRNAs){
    p_temp = mRNA_lncRNA_phyper(m,lnc)
    mRNAs = c(mRNAs,m)
    lncRNAs = c(lncRNAs,lnc)
    phyper_p = c(phyper_p,p_temp)
    cor_coe = c(cor_coe,mRNA_lncRNA_cor$cor[m,lnc])
    cor_p = c(cor_p,mRNA_lncRNA_cor$p[m,lnc])
  }
  cat("have done:",grep(m,diff_mRNAs),m,"\n")
}
mRNA_ce_lncRNA = data.frame(mRNAs,lncRNAs,phyper_p,cor_coe,cor_p,stringsAsFactors = F)
#mRNA_ce_lncRNA = filter(mRNA_ce_lncRNA,p_value<=0.05)
dim(mRNA_ce_lncRNA)
#save.image(file = "step3-ceRNAs.RData")


#############################################################
rm(list = ls())
load("step3-ceRNAs.RData")
mRNA_ce_lncRNA_sig = filter(mRNA_ce_lncRNA, phyper_p<=0.05, cor_coe>=0.4, cor_p<=0.05)
length(unique(mRNA_ce_lncRNA_sig$mRNAs))
length(unique(mRNA_ce_lncRNA_sig$lncRNAs))
####提取ceRNA网络中的节点
diff_miRNAs[13] = "hsa-mir-3199"
ceRNA_miRNA_mRNA = miRNA_mRNA[which(miRNA_mRNA$mRNA %in% mRNA_ce_lncRNA_sig$mRNAs),]
ceRNA_miRNA_mRNA$miRNA = tolower(ceRNA_miRNA_mRNA$miRNA)
ceRNA_miRNA_mRNA = ceRNA_miRNA_mRNA[which(ceRNA_miRNA_mRNA$miRNA %in% diff_miRNAs),]

ceRNA_miRNA_lncRNA = miRNA_lncRNA[which(miRNA_lncRNA$lncRNA %in% mRNA_ce_lncRNA_sig$lncRNAs),]
ceRNA_miRNA_lncRNA$miRNA = tolower(ceRNA_miRNA_lncRNA$miRNA)
ceRNA_miRNA_lncRNA = ceRNA_miRNA_lncRNA[which(ceRNA_miRNA_lncRNA$miRNA %in% diff_miRNAs),]
{
  ceRNA_mRNA_miRNA_lncRNA = inner_join(ceRNA_miRNA_mRNA,ceRNA_miRNA_lncRNA,by="miRNA")
}

index = c()
for(i in 1:dim(ceRNA_mRNA_miRNA_lncRNA)[1]){
  m1 = ceRNA_mRNA_miRNA_lncRNA[i,2]
  l1 = ceRNA_mRNA_miRNA_lncRNA[i,3]
  for(j in 1:dim(mRNA_ce_lncRNA_sig)[1]){
    m2 = mRNA_ce_lncRNA_sig[j,1]
    l2 =mRNA_ce_lncRNA_sig[j,2]
    if(m1==m2 & l1==l2){ index = c(index,i) }
  }
}

ceRNA_network = ceRNA_mRNA_miRNA_lncRNA[index,]
length(unique(ceRNA_network$miRNA))
length(unique(ceRNA_network$mRNA))
length(unique(ceRNA_network$lncRNA))

ceRNA_df = NULL
for(ii in 1:dim(ceRNA_network)[1]){
  ceRNA_df = rbind(ceRNA_df,
                   rbind(c(ceRNA_network[ii,1],ceRNA_network[ii,2],"miRNA","mRNA"),
                         c(ceRNA_network[ii,1],ceRNA_network[ii,3],"miRNA","lncRNA"),
                         c(ceRNA_network[ii,2],ceRNA_network[ii,3],"mRNA","lncRNA")))
}
ceRNA_df = as.data.frame(ceRNA_df)
ceRNA_df = distinct(ceRNA_df)
write.table(ceRNA_df,file = "results/ceRNA_network.txt",sep = "\t",quote = F,row.names = F,col.names = T)

#save.image(file = "step3-ceRNAs22.RData")



sort(table(unlist(ceRNA_df[,1:2])),decreasing = T)

