rm(list = ls())
mRNA1 = c("MCAM")
mRNA = c("DOCK11","MCAM","MYO10","WASF3")
miRNA = c("hsa-mir-141","hsa-mir-143","hsa-mir-140","hsa-mir-139","hsa-mir-145","hsa-mir-5683")

setwd("F:/chenkg/THYM/图表/validation")
geneExp = read.table("GSE79978_gene_RPKM.txt",header = T,row.names = 1,sep = "\t",
                     stringsAsFactors = F,check.names = F)
miRNAExp = read.table("GSE79978_mir_RPKM.txt",header = T,row.names = 1,sep = "\t",
                      stringsAsFactors = F,check.names = F)

normal=6:8
tumor=c(1:5,9:16)

mRNA_test = geneExp[,mRNA]
par(las=1,mar=c(4,5,3,3))
x=c(1:ncol(mRNA_test))
y=c(1:ncol(mRNA_test))
plot(x,y,
     xlim=c(-1,11),ylim=c(min(mRNA_test),max(mRNA_test)*1.1),
     main="GSE79978",xlab="", ylab="mRNA expression",
     pch=21,
     cex.lab=1.5,
     col="white",
     xaxt="n")

library(vioplot)
#对每个基因循环，绘制vioplot，正常用蓝色表示，肿瘤用红色表示
log2FC = c()
pvalue = c()
for(i in 1:ncol(mRNA_test)){
  normalData=mRNA_test[normal,i]
  tumorData=mRNA_test[tumor,i]
  log2FC = c(log2FC,log2(mean(tumorData)/mean(normalData)))
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = '#1A3263')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = '#C82121')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  pvalue = c(pvalue,p)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx*1.05, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.8)
}
text(seq(1,10,3),-0.5,xpd = NA,labels=colnames(mRNA_test),cex = 1,srt = 45,pos=2)
legend("topright", c("Normal","Tumor"), cex = 0.7, fill = c("#1A3263","#C82121"))






miRNA_test = miRNAExp[,miRNA]
par(las=1,mar=c(4,5,3,3))
x=c(1:ncol(miRNA_test))
y=c(1:ncol(miRNA_test))
plot(x,y,
     xlim=c(-1,50),ylim=c(min(miRNA_test),max(miRNA_test)*1.1),
     main="GSE79978",xlab="", ylab="mRNA expression",
     pch=21,
     cex.lab=1.5,
     col="white",
     xaxt="n")

library(vioplot)
#对每个基因循环，绘制vioplot，正常用蓝色表示，肿瘤用红色表示
for(i in 1:ncol(miRNA_test)){
  normalData=miRNA_test[normal,i]
  tumorData=miRNA_test[tumor,i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = '#1A3263')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = '#C82121')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx*1.05, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.8)
}
text(seq(1,50,3),0,xpd = NA,labels=colnames(miRNA_test),cex = 1,srt = 45,pos=2)
legend("top", c("Normal","Tumor"), cex = 0.7, fill = c("#1A3263","#C82121"))





