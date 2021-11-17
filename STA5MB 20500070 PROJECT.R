#loading required packages
library(cluster) 
library(genefilter)
library(gplots)
library(limma)
library(qvalue)
library(RColorBrewer)
library(survival)
library(SingleCellExperiment)
library(scater)
library(irlba)
library(scran)
library(Rtsne)

#loading the data sets
BC_clinical<-readRDS("BC_clinical.rds")
BC_data<-readRDS("BC_data.rds")
colSums(is.na(BC_clinical))
colSums(is.na(BC_data))

#inspect data
head(BC_clinical)
table(BC_clinical$LNstatus)
head(BC_data)

#hierarchial clustering
scaled.BC_data <- t(scale(t(BC_data), center = TRUE))
head(scaled.BC_data)
d <- dist(t(scaled.BC_data))
hc <- hclust(d, method = "complete")

#plot
plot(hc, cex=0.5)
abline(h=220, col = 2)

K <- 2:5
sh <- NULL
for (i in K){
  sh <- c(sh,median(silhouette(cutree(hc,k = i), dist = d)[,3],na.rm = T))
}

#silhouette
plot(K, sh, type = "l", main = "Median silhouette",xlab="num of clusters")

cl = cutree(hc,k=K[which.max(sh)])
table(cl)
#heat map
rv <- rowVars(scaled.BC_data)
idx <- order(-rv)[1:5000]
cols <- colors()[seq(1, length(colors()), len = length(unique(cl)))]
# Inspect colors mapped to columns of BC_data
head(cbind(colnames(BC_data), cols))
# Produce heat-map
suppressWarnings(library(gplots))
heatmap.2(scaled.BC_data[idx, ], labCol = cl, trace = "none",
          ColSideColors = cols[cl],col = redgreen(100),main = 'HEATMAP of CLUSTERS')
suppressWarnings(head(cbind(colnames(BC_data),cols)))
#PCA analysis
par(bg = "grey")
pc <- princomp(scaled.BC_data[idx, ])
plot(pc$load[, 1:2], col = cl)
title("Different sub types of brest cancer data, coloured by clusters")

#Performing the Differential gene expression analysis by comparing the clusters

design <- model.matrix(~as.factor(cl))
DE.object<- lmFit(BC_data, design)
DE.object <- eBayes(DE.object)

qval<- qvalue(DE.object$p.value[,2], fdr.level = 0.05)


#survival
gene.score <- colSums(BC_data[qval$sig, ])
# standardize gene score (to have mean=0, SD=1) 
gene.score <- scale(gene.score)

boxplot(split(gene.score,cl),col=c("yellowgreen","white","grey68"),main="Boxplot
        of gene scores for genes which are DE between clusters")
# Perform Cox regression to estimate HR of gene.score
cox.model <- coxph(Surv(Surv_time, event) ~ gene.score + histgrade + 
                     ERstatus + PRstatus + age + tumor_size_mm + 
                     LNstatus , data = BC_clinical)
summary(cox.model)



