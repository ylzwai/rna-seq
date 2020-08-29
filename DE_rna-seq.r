rm(list = ls())
#BiocManager::install(c("airway","DESeq2","edgeR","limma"))
library(airway)
data("airway")
exprSet <- assay(airway)
group_list <- colData(airway)[,3]

suppressMessages(library(DESeq2)) 
colData <- data.frame(row.names=colnames(exprSet), 
                       group_list=group_list) 
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds)


res <- results(dds, 
               contrast=c("group_list","trt","untrt"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG=as.data.frame(resOrdered)
head(DEG,1)


#DEG里面有很多NA，不能做差异分析。
DEG <- na.omit(DEG)
if(F){
#热图
#install.packages("pheatmap")

library(pheatmap)

need_DEG <- DEG
choose_gene=head(rownames(need_DEG),100) ## 50 maybe better
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename = '_need_DEG_top50_heatmap.png')


#火山图
logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs( log2FoldChange)) )


need_DEG$change <-  as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                   ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',])
)


  library(ggplot2)
g = ggplot(data=need_DEG, 
           aes(x=log2FoldChange, y=-log10(pvalue), 
               color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
print(g)
ggsave(g,filename = 'volcano.png')
}
DEseq_DEG <- DEG

library(edgeR)
d <- DGEList(counts=exprSet,group=factor(group_list))
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d$samples

design <- model.matrix(~0+factor(group_list))
dge <- d
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))

dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)

lrt <- glmLRT(fit,  contrast=c(1,0))
nrDEG=topTags(lrt, n=nrow(exprSet))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
edgeR_DEG <- nrDEG



#limma
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)

dge

dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('untrt-trt'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='untrt-trt', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
##edgeR_DEG代码错误，logFC太大
save(DEG_limma_voom,edgeR_DEG,DEseq_DEG,exprSet,group_list,
     file = "./Rdata/DEG_results.Rdata")


pheatmap()