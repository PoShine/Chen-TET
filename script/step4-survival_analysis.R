rm(list = ls())
library(dplyr)
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library (VennDiagram)
library(data.table)
library(survival)
library(survminer)
library(survivalROC)
setwd("F:/AAA/chenkg/THYM/data")

load("step3-ceRNAs22.RData")

diff_mRNAs_fpkm_t = diff_mRNAs_fpkm[diff_mRNAs,c(colnames(mRNA_carcinoma),colnames(mRNA_typeA),colnames(mRNA_typeB),colnames(mRNA_typeAB))]
dim(diff_mRNAs_fpkm_t)
diff_lncRNAs_fpkm_t = diff_lncRNAs_fpkm[diff_lncRNAs,c(colnames(lncRNA_carcinoma),colnames(lncRNA_typeA),colnames(lncRNA_typeB),colnames(lncRNA_typeAB))]
dim(diff_lncRNAs_fpkm_t)
rownames(diff_miRNAs_RPM)[13] = "hsa-mir-3199"
diff_miRNAs_RPM_t = diff_miRNAs_RPM[diff_miRNAs,colnames(diff_mRNAs_fpkm_t)]
dim(diff_miRNAs_RPM_t)

survival_df = read.table("TCGA-THYM.survival.tsv",header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
survival_df = survival_df[colnames(diff_mRNAs_fpkm_t),c(1,3)]
survival_df$OS.time = survival_df$OS.time/(365/12)

HR=c()
up95=c()
low95=c()
clnm=c()
for(i in 1:dim(ceRNA_network)[1]){
  mirrna_exp = unlist(diff_miRNAs_RPM_t[ceRNA_network[i,1],])
  mrna_exp = unlist(diff_mRNAs_fpkm_t[ceRNA_network[i,2],])
  lncrna_exp = unlist(diff_lncRNAs_fpkm_t[ceRNA_network[i,3],])
  cerna_surv = data.frame(miRNA=mirrna_exp,mRNA=mrna_exp,lncRNA=lncrna_exp,event=survival_df$OS,time=survival_df$OS.time)
  res.cox <- coxph(Surv(time, event) ~ miRNA + mRNA + lncRNA, data =  cerna_surv)
  risk.score = (as.matrix(cerna_surv[,1:3]) %*% matrix(res.cox$coefficients,ncol=1))[,1]
  group <- rep(0, length(risk.score))
  for (j in 1:length(group)) {
    if( risk.score[j]<=median(risk.score)) {
      group[j]<-"low"
    } else{
      group[j]<-"high"
    }
  }
  surv.df <- data.frame(event=survival_df$OS,time=survival_df$OS.time,risk.score,group,stringsAsFactors = F)
  test.survdiff <- survdiff(Surv(time,event)~ group,data=surv.df)
  HR_temp <- (test.survdiff$obs[1]/test.survdiff$exp[1])/(test.survdiff$obs[2]/test.survdiff$exp[2])
  HR_temp <- round(HR_temp,4)
  p.val_temp <- round(1 - pchisq(test.survdiff$chisq, 1),4)
  up95_temp <- exp(log(as.numeric(HR_temp)) + qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  up95_temp <- round(up95_temp,4)
  low95_temp <- exp(log(as.numeric(HR_temp)) - qnorm(0.975)*sqrt(1/test.survdiff$exp[2]+1/test.survdiff$exp[1]))
  low95_temp <- round(low95_temp,4)
  
  surv.fit = survfit(Surv(time, event) ~ group,data = surv.df)
  ptitle = paste(ceRNA_network[i,2],ceRNA_network[i,3],ceRNA_network[i,1],sep = "-")
  ggsurvplot(
    surv.fit,  # survfit object with calculated statistics.
    #surv.median.line = "hv", # Add medians survival
    #risk.table = TRUE,       # show risk table.
    #tables.height = 0.2,
    #tables.theme = theme_cleantable(),
    title = ptitle,
    pval = TRUE,             # show p-value of log-rank test.
    ggtheme = theme_bw(), # customize plot and risk table with a theme.
    #risk.table.y.text.col = T, # colour risk table text annotations.
    #risk.table.y.text = FALSE, # show bars instead of names in text annotations
    palette = c("red","darkblue")
  )
  #ggsave(filename = paste0(ptitle,".pdf"),path = "results/kmcurves/")
  
  
  HR=c(HR,HR_temp)
  up95=c(up95,up95_temp)
  low95=c(low95,low95_temp)
  clnm=c(clnm,ptitle)
  
  #pdf(file=paste0("results/kmcurves/",ptitle,"ROC.pdf"))
  #par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  #roc=survivalROC(Stime=surv.df$time, status=surv.df$event, marker = surv.df$risk.score, 
  #                predict.time = 60, method="KM")
  #plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
  #     xlab="False positive rate", ylab="True positive rate",
  #     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
  #     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  #abline(0,1)
  #dev.off()
}

pre_result=cbind(HR=HR,low95=low95,up95=up95)
rownames(pre_result) = clnm


