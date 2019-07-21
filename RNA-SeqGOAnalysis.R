##R-code used for "Transcriptome analysis of 
##axolotl oropharyngeal explants during taste bud differentiation stages."
################################################################################# 
#---------------------------------------------------------------------------
# this code shows complete analysis used in the paper in following 6 step:
#---------------------------------------------------------------------------
#step1: reading raw counts data and sample information
#step2: pre-processing data
#step3: quality control: MDS plot, PCA plot and dendogram
#step4: correlation analysis
#step5: DEG analysis using Limma
#step6: Go terms for DEGs
#############################################################################
##Author: P.Kohli, Associate Professor of Statistics
##Dept. of Mathematics and Statistics
##Connecticut College
##New London, CT, 06320
##Last Updated on: 7/19/2019
############################################################################

##Load required R packages
library(edgeR)
library(limma)
library(corrplot)
library(locfit)
#######################################################################
## Step1: Reading Data ##
#######################################################################
##reading raw counts and sample information
data <- read.table("Second2017_counts.counts.matrix", header=T, row.names=1, com='')
sampleinfo <- read.table("sampleinfo.txt",header = TRUE)
names(data)
##naming stages and replicates
names(data) <- c("37K","37L","37M","38A","38B","38C","39D","39E","39F","39G","39N","39O","41H","41J","41P","41Q","41R")
##removing 43S and T for the quality control analysis
data <- data[,-c(18,19)]
#remove corresponding stage from sampleinfo
sampleinfo <- sampleinfo[-c(18,19),]

#######################################################################
## Step2: Pre-processing Data ##
#######################################################################
#counts per million
CPM <- cpm(data)
#keep genes that have at least two TRUEs in each row of thresh
thresh <- CPM>1
keep <- rowSums(thresh)>=2
counts.keep <- data[keep,]
#Creates a DGEList object from a table of counts, 
##library size information is stored in the sample lot
y <- DGEList(counts.keep)
##plot lib sizes as bar plots to see any major discrepencies
barplot(y$samples$lib.size,names=colnames(y),las=2,xlab="Sample",ylab="Library Size",cex.axis = 0.7,cex.names = 0.7,yaxt="n")
axis(2,at=c(0,100000,200000,300000,400000,500000,600000,700000),labels=c("0","100,000","200,000","300,000","400,000","500,000","600,000","700,000"),cex.axis=0.7)
title("Barplot of library sizes")
##log of counts data is used 
## cpm function can return log2cpm which are corrected for 
##library sizes and adds a small offset to avoid taking log of zero
y <- calcNormFactors(y)
#get log2 counts per million
logcounts <- cpm(y,log=TRUE)
barplot(logcounts,names=colnames(y),las=2,xlab="Sample",ylab="Library Size",cex.axis = 0.7,cex.names = 0.7)
normalized.counts <- cpm(y)
transposed <- t(logcounts)

#######################################################################
## Step3: Quality Control ##
#######################################################################
##1. MDSplot 
samplR <- c(rep("37/38",6),rep(39,6),rep(41,5))
levels(as.factor(samplR))

col.cell<-c(1,2,3)[as.factor(samplR)]
data.frame(as.factor(samplR),col.cell)
plotMDS(y,col=col.cell,gene.selection = "pairwise",cex=1)
legend("top",horiz=T,legend=levels(as.factor(samplR)),fill=c(1,2,3,4,5),cex=1)
title("(a)")

##2. PCA Plot
pca <- prcomp(t(logcounts), scale = T)
plot(pca)
summary(pca)
loadings <- pca$rotation
scores <- pca$x

plot(scores[,1], scores[,2], main="(b)",
     xlab=paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained"), 
     ## indicate the % variance explained by PC1
     ylab=paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained"), 
     ## indicate the % variance explained by PC2
     pch=16,col=col.cell,cex=1.5,
     ## points shape according to treatment
     xlim=c(min(scores[,1]), max(scores[,1])) ,
     ylim=c(min(scores[,2]), max(scores[,2])+(max(scores[,2])-min(scores[,2]))/4) 
     ## Let a bit of room on top of the plot for legend
)
data2 <- data
names(data2)[1:6]<-c("37K", "37L","37M","38A","38B","38C")
text(scores[,1], scores[,2] + 5, names(data2),cex=0.8)
legend("top", horiz=T, legend=levels(as.factor(samplR)),fill=c(1,2,3,4,5),cex=0.8)

##3. Dendogram
distance <- dist(transposed)
clusters <- hclust(distance)
plot(clusters,main="(c)",ylab="Euclidean Distance",xlab=NA,sub=NA)

##All quality control plots together
par(mfrow=c(1,3))
##MDS
samplR<-c(rep("37/38",6),rep(39,6),rep(41,5))
levels(as.factor(samplR))
col.cell<-c(1,2,3)[as.factor(samplR)]
data.frame(as.factor(samplR),col.cell)
plotMDS(y,col=col.cell,gene.selection = "pairwise",cex=1.2)
legend("top",horiz=T,legend=levels(as.factor(samplR)),fill=c(1,2,3,4,5),cex=1)
title("(a)")
##PCA
plot(scores[,1], scores[,2], main="(b)",
     xlab=paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained"), 
     ## indicate the % variance explained by PC1
     ylab=paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained"), 
     ## indicate the % variance explained by PC2
     pch=16,col=col.cell,cex=1.5,
     ## points shape according to treatment
     xlim=c(min(scores[,1]), max(scores[,1])) ,
     ylim=c(min(scores[,2]), max(scores[,2])+(max(scores[,2])-min(scores[,2]))/4) 
     ## Let a bit of room on top of the plot for legend
)
text(scores[,1], scores[,2] + 5, names(data2),cex=1.2)
legend("top", horiz=T, legend=levels(as.factor(samplR)),fill=c(1,2,3,4,5),cex=1)
##Dendogram
plot(clusters,main="(c)",ylab="Euclidean Distance",xlab=NA,sub=NA)

#######################################################################
## Step4: Pearson's Correlations ##
#######################################################################

##Correlation matrices
logcounts <- as.data.frame(logcounts)

##correlaion matrices
stage1.cor <- cor(logcounts[,1:6])
stage2.cor <- cor(logcounts[,7:12])
stage3.cor <- cor(logcounts[,13:17])

par(mfrow=c(1,3))
corrplot(stage1.cor,method="number",tl.col="red",cl.lim=c(0,1),title="(a)",type="upper",diag=FALSE,mar=c(0,0,1,0))
corrplot(stage2.cor,method="number",tl.col="red",cl.lim=c(0,1),title="(b)",type="upper",diag=FALSE,mar=c(0,0,1,0))
corrplot(stage3.cor,method="number",tl.col="red",cl.lim=c(0,1),title="(c)",type="upper",diag=FALSE,mar=c(0,0,1,0))


#######################################################################
## Step5: DEG Analysis using Limma for selected replicates comparisons##
#######################################################################
##Read matrix with counts data
data <- read.table("Second2017_counts.counts.matrix", header=T, row.names=1, com='')
names(data)
##select replicates to be compared
#For instance: 39FO versus 41QR
col_ordering <- c(9,12,16,17)
###select subset of data from stages
names(data)[col_ordering]
rnaseqMatrix <- data[,col_ordering]
##rounding to 0 decimals
rnaseqMatrix <- round(rnaseqMatrix)
##filtering low counts
rnaseqMatrix <- rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
##setting contrast and design 
conditions <- factor(c(rep("39", 2), rep("41", 2)))

##Limma: design matrix
design <- model.matrix(~0+conditions)
##name of stages
colnames(design) <- c("Stage39","Stage41")
##make contrast matrix
contrast.matrix <- makeContrasts(Stage41-Stage39,levels=design)
## TMM normalize data: to align columns of a count matrix
lib_sizes <- colSums(rnaseqMatrix)
tmm_norm_factors <- calcNormFactors(rnaseqMatrix, method='TMM')
##MAkes a DGE List from counts and group indicators
x <- DGEList(counts=rnaseqMatrix,group = conditions)
# voom transformation for linear modeling of the counts data
y <- voom(x, design, lib.size=lib_sizes*tmm_norm_factors, plot=T)
##Empirical Bayes for DE Analysis
fit <- eBayes(contrasts.fit(lmFit(y,design),contrast.matrix))
## extracts top-ranked genes from above result
tTags <- topTable(fit,coef=1,number=Inf)
nrow(tTags[tTags$adj.P.Val < 0.05, ])
##select DEGs using significant adjusted P-values 
tTags2 <- tTags[tTags$adj.P.Val < 0.05, ]
##check count of positive logFC: which tells
sum(tTags2$logFC >0)
##check count of negative logFC: which tells
sum(tTags2$logFC <0)
##conduct multiple testing across genes and contrasts
results <- decideTests(fit)
##VennDiagram
vennCounts(results)

##Volcano Plot
##counts per million for x
c <- cpm(x)
##row mean
m <- apply(c, 1, mean)
##log FC
tTags$logFC <- tTags$logFC  
#setting results
tTags2 <- cbind(tTags, logCPM=log2(m[rownames(tTags)]))
DE_matrix <- data.frame(sampleA="39FO", sampleB="41QR", logFC=tTags$logFC, logCPM=tTags2$logCPM, PValue=tTags$'P.Value', FDR=tTags$'adj.P.Val')
rownames(DE_matrix) <- rownames(tTags)
DE_matrix <- data.frame(sampleA="39FO", sampleB="41QR", logFC=tTags$logFC, logCPM=tTags2$logCPM, PValue=tTags$'P.Value', FDR=tTags$'adj.P.Val')
rownames(DE_matrix) <- rownames(tTags)

plot_Volcano2 = function(data, xlab="log2(FoldChange)", ylab="-log10(FDR)", title="", pch=20,legendtext) {
  par(mar = c(5, 4, 2, 0.2)) 
  subset.data <- data[data$adj.P.Val<0.05,]
  plot(data$logFC, -1*log10(data$adj.P.Val), col=ifelse(data$adj.P.Val<0.05 & data$logFC>0,"red",ifelse(data$adj.P.Val<0.05 & data$logFC<0,"orange","gray")), xlab=xlab, ylab=ylab, main=title, pch=pch,cex.lab=0.7)
  legend("topright",legend=legendtext,horiz=FALSE,col=c("red","orange","gray"),pch=20,cex=0.9)
}

plot_Volcano2(tTags,title="(i):39FO vs 41QR",legendtext=c("Upregulated in Stage 41","Upregulated in Stage 39","No Difference"))

##counts for DEGs at different threshold levels of logFC
DEG <- tTags[tTags$adj.P.Val < 0.05, ]
##all up in 41
nrow(DEG[DEG$logFC>0,])
#all up in 39
nrow(DEG[DEG$logFC<0,])
##DEGs
nrow(DEG[abs(DEG$logFC)>1,])
nrow(DEG[abs(DEG$logFC)>2,])
nrow(DEG[abs(DEG$logFC)>3,])
nrow(DEG[abs(DEG$logFC)>4,])
nrow(DEG[abs(DEG$logFC)>5,])
##up in 41 at different levels
nrow(DEG[DEG$logFC>1,])
nrow(DEG[DEG$logFC>2,])
nrow(DEG[DEG$logFC>3,])
nrow(DEG[DEG$logFC>4,])
nrow(DEG[DEG$logFC>5,])
##up in 39 at different levels
nrow(DEG[abs(DEG$logFC)>1,])-nrow(DEG[DEG$logFC>1,])
nrow(DEG[abs(DEG$logFC)>2,])-nrow(DEG[DEG$logFC>2,])
nrow(DEG[abs(DEG$logFC)>3,])-nrow(DEG[DEG$logFC>3,])
nrow(DEG[abs(DEG$logFC)>4,])-nrow(DEG[DEG$logFC>4,])
nrow(DEG[abs(DEG$logFC)>5,])-nrow(DEG[DEG$logFC>5,])


#######################################################################
## Step6: GO Terms identification##
#######################################################################
##GO Annotation: matching
##read blast2go created annotation table
annotations <- read.csv("blast2go_annotated_table.csv",header=TRUE)
##check for duplication
sum(duplicated(annotations[,2]))
seqnames <- as.character(annotations[,2])
step1 <- noquote(substr(seqnames,start=1,stop=21))
annotations[,2] <- step1
sum(duplicated(annotations[,2]))
results <- DEG
sum(duplicated(results[,1]))

write.csv(results,"LimmaDEG-39FN41HP.csv")
DEG <- read.csv("LimmaDEG-39FN41HP.csv")
names(DEG)[1]<-"SeqName"
final.data <- merge(annotations,DEG,by.x = "SeqName")
check <- annotations
check2 <- check[with(check, SeqName %in% results$SeqName), ]
##removing duplicates
annotations.new <-annotations[!duplicated(annotations[,2]),] 
final.data.new <- merge(annotations.new,DEG,by.x = "SeqName")
##saving final DEG results with GO terms
write.csv(final.data.new,"ResultsComparisonLimmaDEG-39FNvs41HP.csv")

