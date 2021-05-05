rm(list = ls())
library(dplyr)
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library (VennDiagram)
setwd("F:/chenkg/THYM/data")
##读取表型、临床数据
phenotype = read.table("TCGA-THYM.GDC_phenotype.tsv",header = T,row.names = 1,quote = "",
                       sep = "\t",stringsAsFactors = F,check.names = F)
phenotype = phenotype[grep("-01A",rownames(phenotype)),]
clinical = read.csv("TCGA-THYM-clinical.csv",header = T,row.names = 1,stringsAsFactors = F)
clinical = clinical[-grep("_1",rownames(clinical)),]
##确定每个样本的疾病亚型
sample_type = phenotype$primary_diagnosis.diagnoses
sample_type = sub("Thymoma, ","",sample_type)
sample_type = sub(", (malignant|NOS)","",sample_type)
sample_type[grep("type B",sample_type)] = "type B"
table(sample_type)
names(sample_type) = rownames(phenotype)
{
  sample_carcinoma = names(sample_type[which(sample_type=="Thymic carcinoma")])
  sample_typeA = names(sample_type[which(sample_type=="type A")])
  sample_typeB = names(sample_type[which(sample_type=="type B")])
  sample_typeAB = names(sample_type[which(sample_type=="type AB")])
}

#导入表达数据，识别DEGs
load("mRNAandLncRNAandmiRNAexp.RData")
##保留至少在10%样本中有表达的基因
ratio = function(mat){
  return(apply(mat,1,
               function(x){
                 return(length(which(x>0))/length(x))
               }))
}
##mRNA
mRNA_ratio = ratio(mRNA_expr_fpkm)
mRNA_exp_counts = mRNA_exp_counts[which(mRNA_ratio>0.1),]
mRNA_expr_fpkm = mRNA_expr_fpkm[which(mRNA_ratio>0.1),]
##lncRNA
lncRNA_ratio = ratio(lncRNA_exp_fpkm)
lncRNA_exp_counts = lncRNA_exp_counts[which(lncRNA_ratio>0.1),]
lncRNA_exp_fpkm = lncRNA_exp_fpkm[which(lncRNA_ratio>0.1),]
##miRNA
miRNA_ratio = ratio(miRNA_exp_RPM)
miRNA_exp_count = miRNA_exp_count[which(miRNA_ratio>0.1),]
miRNA_exp_RPM = miRNA_exp_RPM[which(miRNA_ratio>0.1),]

###########mRNA
if(T){
  #mRNA 
  dim(mRNA_exp_counts)
  mRNA_carcinoma = mRNA_exp_counts[,intersect(colnames(mRNA_exp_counts),sample_carcinoma)]
  dim(mRNA_carcinoma)
  mRNA_typeA = mRNA_exp_counts[,intersect(colnames(mRNA_exp_counts),sample_typeA)]
  dim(mRNA_typeA)
  mRNA_typeB = mRNA_exp_counts[,intersect(colnames(mRNA_exp_counts),sample_typeB)]
  dim(mRNA_typeB)
  mRNA_typeAB = mRNA_exp_counts[,intersect(colnames(mRNA_exp_counts),sample_typeAB)]
  dim(mRNA_typeAB)
  mRNA_normal = mRNA_exp_counts[,grep("-11A",colnames(mRNA_exp_counts))]
  dim(mRNA_normal)
  
  ##DEGs
  mRNA_DEGs = function(stype,alpha=0.25){
    mRNA_normal_stype = cbind(mRNA_normal,eval(parse(text = paste0("mRNA_",stype))))
    mRNA_normal_stype = round(mRNA_normal_stype)
    #sample info
    condition <- factor(c(rep("normal",2),rep(stype,dim(eval(parse(text = paste0("mRNA_",stype))))[2])), levels = c("normal",stype))  #condition是因子
    colData <- data.frame(row.names=colnames(mRNA_normal_stype), condition)
    dds <- DESeqDataSetFromMatrix(mRNA_normal_stype, colData, design= ~ condition)
    dds <- DESeq(dds)
    #查看stype versus normal的总体结果，并根据fdr进行重新排序。
    #利用summary命令统计显示一共多少个genes上调和下调（FDR）
    res = results(dds, contrast=c("condition", stype, "normal"),alpha = alpha,pAdjustMethod = "fdr")
    res = res[order(res$padj),]
    head(res)
    summary(res)
    #获取padj 小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因
    res.df <- na.omit(as.data.frame(res))
    diff_mRNA_stype <- subset(res.df,padj < alpha & (log2FoldChange > 1 | log2FoldChange < -1))
    return(diff_mRNA_stype)
  }
  
  #DEGs carcinoma vs. normal
  diff_mRNA_carcinoma = mRNA_DEGs("carcinoma")
  diff_mRNA_typeA = mRNA_DEGs("typeA")
  diff_mRNA_typeB = mRNA_DEGs("typeB")
  diff_mRNA_typeAB = mRNA_DEGs("typeAB")
  
  demRNAs = intersect(intersect(intersect(rownames(diff_mRNA_carcinoma),rownames(diff_mRNA_typeA)),rownames(diff_mRNA_typeB)),rownames(diff_mRNA_typeAB))
  upmRNA1 = demRNAs[which(diff_mRNA_carcinoma[demRNAs,]$log2FoldChange > 0)]
  upmRNA2 = demRNAs[which(diff_mRNA_typeA[demRNAs,]$log2FoldChange > 0)]
  upmRNA3 = demRNAs[which(diff_mRNA_typeB[demRNAs,]$log2FoldChange > 0)]
  upmRNA4 = demRNAs[which(diff_mRNA_typeAB[demRNAs,]$log2FoldChange > 0)]
  upmRNA = intersect(intersect(upmRNA1,upmRNA2),intersect(upmRNA3,upmRNA4))
  
  downmRNA1 = demRNAs[which(diff_mRNA_carcinoma[demRNAs,]$log2FoldChange < 0)]
  downmRNA2 = demRNAs[which(diff_mRNA_typeA[demRNAs,]$log2FoldChange < 0)]
  downmRNA3 = demRNAs[which(diff_mRNA_typeB[demRNAs,]$log2FoldChange < 0)]
  downmRNA4 = demRNAs[which(diff_mRNA_typeAB[demRNAs,]$log2FoldChange < 0)]
  downmRNA = intersect(intersect(downmRNA1,downmRNA2),intersect(downmRNA3,downmRNA4))
  
  delmRNA = setdiff(demRNAs,c(upmRNA,downmRNA))
  diff_mRNA_carcinoma = diff_mRNA_carcinoma[-which(rownames(diff_mRNA_carcinoma) %in% delmRNA),]
  diff_mRNA_typeA = diff_mRNA_typeA[-which(rownames(diff_mRNA_typeA) %in% delmRNA),]
  diff_mRNA_typeB = diff_mRNA_typeB[-which(rownames(diff_mRNA_typeB) %in% delmRNA),]
  diff_mRNA_typeAB = diff_mRNA_typeAB[-which(rownames(diff_mRNA_typeAB) %in% delmRNA),]
  
  ##使用venn.diagram功能绘图
  venn.diagram(x= list('Thymic carcinoma\n_vs_Normal' = rownames(diff_mRNA_carcinoma),"typeA_vs_Normal" = rownames(diff_mRNA_typeA),"typeB_vs_Normal" = rownames(diff_mRNA_typeB),"typeAB_vs_Normal" = rownames(diff_mRNA_typeAB)), 
               filename = "plots/mRNA_DEGs.png",height = 1200, width = 1200,resolution =300, imagetype="png", col="transparent",
               fill=c("cornflowerblue","gray","orangered2","darkorchid1"),alpha = 0.4, cex=0.45, cat.cex=0.5)

}

###########lncRNA
if(T){
  #lncRNA
  dim(lncRNA_exp_counts)
  lncRNA_carcinoma = lncRNA_exp_counts[,intersect(colnames(lncRNA_exp_counts),sample_carcinoma)]
  dim(lncRNA_carcinoma)
  lncRNA_typeA = lncRNA_exp_counts[,intersect(colnames(lncRNA_exp_counts),sample_typeA)]
  dim(lncRNA_typeA)
  lncRNA_typeB = lncRNA_exp_counts[,intersect(colnames(lncRNA_exp_counts),sample_typeB)]
  dim(lncRNA_typeB)
  lncRNA_typeAB = lncRNA_exp_counts[,intersect(colnames(lncRNA_exp_counts),sample_typeAB)]
  dim(lncRNA_typeAB)
  lncRNA_normal = lncRNA_exp_counts[,grep("-11A",colnames(lncRNA_exp_counts))]
  dim(lncRNA_normal)
  
  ##DEGs
  lncRNA_DEGs = function(stype,alpha=0.25){
    lncRNA_normal_stype = cbind(lncRNA_normal,eval(parse(text = paste0("lncRNA_",stype))))
    lncRNA_normal_stype = round(lncRNA_normal_stype)
    #sample info
    condition <- factor(c(rep("normal",2),rep(stype,dim(eval(parse(text = paste0("lncRNA_",stype))))[2])), levels = c("normal",stype))  #condition是因子
    colData <- data.frame(row.names=colnames(lncRNA_normal_stype), condition)
    dds <- DESeqDataSetFromMatrix(lncRNA_normal_stype, colData, design= ~ condition)
    dds <- DESeq(dds)
    #查看stype versus normal的总体结果，并根据fdr进行重新排序。
    #利用summary命令统计显示一共多少个genes上调和下调（FDR）
    res = results(dds, contrast=c("condition", stype, "normal"),alpha = alpha,pAdjustMethod = "fdr")
    res = res[order(res$padj),]
    head(res)
    summary(res)
    #获取padj 小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因
    res.df <- na.omit(as.data.frame(res))
    diff_lncRNA_stype <- subset(res.df,padj < alpha & (log2FoldChange > 1 | log2FoldChange < -1))
    return(diff_lncRNA_stype)
  }
  
  #DEGs carcinoma vs. normal
  diff_lncRNA_carcinoma = lncRNA_DEGs("carcinoma")
  diff_lncRNA_typeA = lncRNA_DEGs("typeA")
  diff_lncRNA_typeB = lncRNA_DEGs("typeB")
  diff_lncRNA_typeAB = lncRNA_DEGs("typeAB")
  
  delncRNAs = intersect(intersect(intersect(rownames(diff_lncRNA_carcinoma),rownames(diff_lncRNA_typeA)),rownames(diff_lncRNA_typeB)),rownames(diff_lncRNA_typeAB))
  uplncRNA1 = delncRNAs[which(diff_lncRNA_carcinoma[delncRNAs,]$log2FoldChange > 0)]
  uplncRNA2 = delncRNAs[which(diff_lncRNA_typeA[delncRNAs,]$log2FoldChange > 0)]
  uplncRNA3 = delncRNAs[which(diff_lncRNA_typeB[delncRNAs,]$log2FoldChange > 0)]
  uplncRNA4 = delncRNAs[which(diff_lncRNA_typeAB[delncRNAs,]$log2FoldChange > 0)]
  uplncRNA = intersect(intersect(uplncRNA1,uplncRNA2),intersect(uplncRNA3,uplncRNA4))
  
  downlncRNA1 = delncRNAs[which(diff_lncRNA_carcinoma[delncRNAs,]$log2FoldChange < 0)]
  downlncRNA2 = delncRNAs[which(diff_lncRNA_typeA[delncRNAs,]$log2FoldChange < 0)]
  downlncRNA3 = delncRNAs[which(diff_lncRNA_typeB[delncRNAs,]$log2FoldChange < 0)]
  downlncRNA4 = delncRNAs[which(diff_lncRNA_typeAB[delncRNAs,]$log2FoldChange < 0)]
  downlncRNA = intersect(intersect(downlncRNA1,downlncRNA2),intersect(downlncRNA3,downlncRNA4))
  
  dellncRNA = setdiff(delncRNAs,c(uplncRNA,downlncRNA))
  diff_lncRNA_carcinoma = diff_lncRNA_carcinoma[-which(rownames(diff_lncRNA_carcinoma) %in% dellncRNA),]
  diff_lncRNA_typeA = diff_lncRNA_typeA[-which(rownames(diff_lncRNA_typeA) %in% dellncRNA),]
  diff_lncRNA_typeB = diff_lncRNA_typeB[-which(rownames(diff_lncRNA_typeB) %in% dellncRNA),]
  diff_lncRNA_typeAB = diff_lncRNA_typeAB[-which(rownames(diff_lncRNA_typeAB) %in% dellncRNA),]
  
  ##使用venn.diagram功能绘图
  venn.diagram(x= list('Thymic carcinoma\n_vs_Normal' = rownames(diff_lncRNA_carcinoma),"typeA_vs_Normal" = rownames(diff_lncRNA_typeA),"typeB_vs_Normal" = rownames(diff_lncRNA_typeB),"typeAB_vs_Normal" = rownames(diff_lncRNA_typeAB)), 
               filename = "plots/lncRNA_DEGs.png",height = 1200, width = 1200,resolution =300, imagetype="png", col="transparent",
               fill=c("cornflowerblue","gray","orangered2","darkorchid1"),alpha = 0.4, cex=0.45, cat.cex=0.5)

}

###########miRNA
if(T){
  #miRNA
  dim(miRNA_exp_count)
  miRNA_carcinoma = miRNA_exp_count[,intersect(colnames(miRNA_exp_count),colnames(mRNA_carcinoma))]
  dim(miRNA_carcinoma)
  miRNA_typeA = miRNA_exp_count[,intersect(colnames(miRNA_exp_count),colnames(mRNA_typeA))]
  dim(miRNA_typeA)
  miRNA_typeB = miRNA_exp_count[,intersect(colnames(miRNA_exp_count),colnames(mRNA_typeB))]
  dim(miRNA_typeB)
  miRNA_typeAB = miRNA_exp_count[,intersect(colnames(miRNA_exp_count),colnames(mRNA_typeAB))]
  dim(miRNA_typeAB)
  miRNA_normal = miRNA_exp_count[,grep("-11A",colnames(miRNA_exp_count))]
  dim(miRNA_normal)
  
  ##DEGs
  miRNA_DEGs = function(stype,alpha=0.25){
    miRNA_normal_stype = cbind(miRNA_normal,eval(parse(text = paste0("miRNA_",stype))))
    miRNA_normal_stype = round(miRNA_normal_stype)
    #sample info
    condition <- factor(c(rep("normal",2),rep(stype,dim(eval(parse(text = paste0("miRNA_",stype))))[2])), levels = c("normal",stype))  #condition是因子
    colData <- data.frame(row.names=colnames(miRNA_normal_stype), condition)
    dds <- DESeqDataSetFromMatrix(miRNA_normal_stype, colData, design= ~ condition)
    dds <- DESeq(dds)
    #查看stype versus normal的总体结果，并根据fdr进行重新排序。
    #利用summary命令统计显示一共多少个genes上调和下调（FDR）
    res = results(dds, contrast=c("condition", stype, "normal"),alpha = alpha,pAdjustMethod = "fdr")
    res = res[order(res$padj),]
    head(res)
    summary(res)
    #获取padj 小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因
    res.df <- na.omit(as.data.frame(res))
    diff_miRNA_stype <- subset(res.df,padj < alpha & (log2FoldChange > 1 | log2FoldChange < -1))
    return(diff_miRNA_stype)
  }
  
  #DEGs carcinoma vs. normal
  diff_miRNA_carcinoma = miRNA_DEGs("carcinoma")
  diff_miRNA_typeA = miRNA_DEGs("typeA")
  diff_miRNA_typeB = miRNA_DEGs("typeB")
  diff_miRNA_typeAB = miRNA_DEGs("typeAB")
  
  demiRNAs = intersect(intersect(intersect(rownames(diff_miRNA_carcinoma),rownames(diff_miRNA_typeA)),rownames(diff_miRNA_typeB)),rownames(diff_miRNA_typeAB))
  upmiRNA1 = demiRNAs[which(diff_miRNA_carcinoma[demiRNAs,]$log2FoldChange > 0)]
  upmiRNA2 = demiRNAs[which(diff_miRNA_typeA[demiRNAs,]$log2FoldChange > 0)]
  upmiRNA3 = demiRNAs[which(diff_miRNA_typeB[demiRNAs,]$log2FoldChange > 0)]
  upmiRNA4 = demiRNAs[which(diff_miRNA_typeAB[demiRNAs,]$log2FoldChange > 0)]
  upmiRNA = intersect(intersect(upmiRNA1,upmiRNA2),intersect(upmiRNA3,upmiRNA4))
  
  ##使用venn.diagram功能绘图
  venn.diagram(x= list('Thymic carcinoma\n_vs_Normal' = rownames(diff_miRNA_carcinoma),"typeA_vs_Normal" = rownames(diff_miRNA_typeA),"typeB_vs_Normal" = rownames(diff_miRNA_typeB),"typeAB_vs_Normal" = rownames(diff_miRNA_typeAB)), 
               filename = "plots/miRNA_DEGs.png",height = 1200, width = 1200,resolution =300, imagetype="png", col="transparent",
               fill=c("cornflowerblue","gray","orangered2","darkorchid1"),alpha = 0.4, cex=0.45, cat.cex=0.5)
  
  
}

#save.image(file = "step2-DEGs.RData")
###############差异基因热图
library(pheatmap)

####mRNA pheatmap
if(T){
  #mRNA, mRNA_expr_fpkm
  diff_mRNAs = intersect(intersect(rownames(diff_mRNA_carcinoma),rownames(diff_mRNA_typeA)),
                         intersect(rownames(diff_mRNA_typeB),rownames(diff_mRNA_typeAB)))
  diff_mRNAs_fpkm = mRNA_expr_fpkm[diff_mRNAs,c(colnames(mRNA_carcinoma),colnames(mRNA_typeA),colnames(mRNA_typeB),colnames(mRNA_typeAB),colnames(mRNA_normal))]
  dim(diff_mRNAs_fpkm)
  
  #col color bar
  annotation_col = data.frame(
    Type = factor(rep(c("Tumor","Normal"), c(119,2)))
  )
  rownames(annotation_col) <- colnames(diff_mRNAs_fpkm)
  # Specify colors
  ann_colors = list(
    Type = c(Tumor = "#f0134d", Normal = "#315b96")
  )
  
  pheatmap(diff_mRNAs_fpkm, cluster_rows = T, cluster_cols = T, show_colnames = F, show_rownames = F,
           scale = "row", annotation_col = annotation_col, annotation_colors = ann_colors, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

####lncRNA pheatmap
if(T){
  #lncRNA, lncRNA_exp_fpkm
  diff_lncRNAs = intersect(intersect(rownames(diff_lncRNA_carcinoma),rownames(diff_lncRNA_typeA)),
                         intersect(rownames(diff_lncRNA_typeB),rownames(diff_lncRNA_typeAB)))
  diff_lncRNAs_fpkm = lncRNA_exp_fpkm[diff_lncRNAs,c(colnames(lncRNA_carcinoma),colnames(lncRNA_typeA),colnames(lncRNA_typeB),colnames(lncRNA_typeAB),colnames(lncRNA_normal))]
  dim(diff_lncRNAs_fpkm)
  
  #col color bar
  annotation_col = data.frame(
    Type = factor(rep(c("Tumor","Normal"), c(119,2)))
  )
  rownames(annotation_col) <- colnames(diff_lncRNAs_fpkm)
  # Specify colors
  ann_colors = list(
    Type = c(Tumor = "#f0134d", Normal = "#315b96")
  )
  
  pheatmap(diff_lncRNAs_fpkm, cluster_rows = T, cluster_cols = T, show_colnames = F, show_rownames = F,
           scale = "row", annotation_col = annotation_col, annotation_colors = ann_colors, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

####miRNA pheatmap
if(T){
  #miRNA, miRNA_exp_RPM
  diff_miRNAs = intersect(intersect(rownames(diff_miRNA_carcinoma),rownames(diff_miRNA_typeA)),
                         intersect(rownames(diff_miRNA_typeB),rownames(diff_miRNA_typeAB)))
  diff_miRNAs_RPM = miRNA_exp_RPM[diff_miRNAs,c(colnames(miRNA_carcinoma),colnames(miRNA_typeA),colnames(miRNA_typeB),colnames(miRNA_typeAB),colnames(miRNA_normal))]
  dim(diff_miRNAs_RPM)
  
  #col color bar
  annotation_col = data.frame(
    Type = factor(rep(c("Tumor","Normal"), c(119,2)))
  )
  rownames(annotation_col) <- colnames(diff_miRNAs_RPM)
  # Specify colors
  ann_colors = list(
    Type = c(Tumor = "#f0134d", Normal = "#315b96")
  )
  
  pheatmap(diff_miRNAs_RPM, cluster_rows = T, cluster_cols = T, show_colnames = F, show_rownames = F,
           scale = "row", annotation_col = annotation_col, annotation_colors = ann_colors, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

########Functional enrichment of mRNA DEGs
diff_mRNAs
library(clusterProfiler)
#symbol转为ENTREZID
gene.df <- bitr(diff_mRNAs, fromType="SYMBOL",
                toType="ENTREZID", 
                OrgDb = "org.Hs.eg.db")
go <- enrichGO(gene = as.numeric(gene.df$ENTREZID), OrgDb = "org.Hs.eg.db", ont="all")
dim(go)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")

kegg <- enrichKEGG(gene = as.numeric(gene.df$ENTREZID), organism = "hsa", keyType = "kegg", 
                   pAdjustMethod = "none", use_internal_data = T)
dim(kegg)
dotplot(kegg, showCategory=30)

save.image(file = "step2-DEGs.RData")
