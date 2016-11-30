######AMEYA KULKARNI#############
######MILES RNASEQ ANALYSIS######
setwd("Desktop/MILES/MILES-surface/Files/")

## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocLite()

#MILES Experimental design
#MUSCLE

filenames <- read.table(file="muscle-files.txt") #muscle-files.txt is a list of IDs submitted to the sequencing core
filenames <- as.character(filenames$V1)
filenames <- strsplit(filenames, split = "[.]")
filenames <- lapply(filenames, function(x){x[1]})
filenames<-filenames[-c(57,47,48)]
unlist(filenames)

subjects <- strsplit(as.character(filenames), split="[-]")
subjects <- lapply(subjects, function(x){x[2]})
unlist(subjects)

treatment <- strsplit(as.character(filenames), split="[-]")
treatment <- lapply(treatment, function(x){x[4]})
unlist(treatment)

visit <- strsplit(as.character(filenames), split="[-]")
visit <- lapply(visit, function(x){x[3]})
unlist(visit)

targets <- cbind(filenames, subjects, treatment, visit)
targets<- data.frame(targets)

library(edgeR)

##Read mus-matrix file created from rsem
muscle <- read.table(file="mus-matrix.txt", header=TRUE) #mus-matrix.txt is the RSEM-generated expected counts matrix

#Round off expected counts to the nearest lower integer
muscle <- round(muscle, digits=0)

#Data matrix for expected counts without columns containing subjects with one biopsy
muscle <- muscle[-c(57,47,48)] #Use this counts matrix for all further analysis

#Generating groups for edgeR
treatment <- factor(as.character(targets$treatment))
subjects <- factor(as.character(targets$subjects))
visit <- factor(as.character(targets$visit))

#Generating a DGEList specific for edgeR analysis in R
genelist <- DGEList(counts=muscle, group=treatment)
genelist <- calcNormFactors(genelist)
cpmlist <- cpm(genelist$counts)
logcpm <- log(cpm(genelist$counts))

#Filtering
keep <- rowSums((cpm(genelist))>1) >=2
genelist <- genelist[keep, , keep.lib.sizes=FALSE]


#Design restructuring

design <- model.matrix(~treatment+subjects+visit)
y <- estimateGLMCommonDisp(genelist,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)


fit <- glmFit(y,design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

plotBCV(y)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

#Gene names
genenames <- strsplit(as.character(rownames((topTags(lrt, n= Inf)))), split="[_]")
genenames <- lapply(genenames, function(x){x[2]})

degtsv <- topTags(lrt, n=Inf)
genenames <- data.frame(unlist(genenames))
degtsv <- cbind(genenames, degtsv)
write.csv (degtsv, file="degtsv.csv")

#Summary tables
libsize<-(colSums(muscle))/(1e06) #library sizes in millions of reads

#volcano Plots
library(calibrate)
allgenes <- read.delim(file = "allgenes.txt", header = TRUE)
genesym <- allgenes$Gene.symbol
#for FC>2 or FC <0.5
with(allgenes, plot(logFC, -log10(PValue), pch=20, main="DEG in muscle after Metformin treatment", xlim=c(-4,5)))
with(subset(allgenes, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(allgenes, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(allgenes, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
with(subset(allgenes, FDR<.05 & abs(logFC)>1), textxy(logFC, -log10(PValue), cex=.6))

#for FC>1.5 or FC <(1/1.5)
with(allgenes, plot(logFC, -log10(PValue), pch=20, main="DEG in muscle after Metformin treatment", xlim=c(-4,5)))
with(subset(allgenes, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(allgenes, abs(logFC)>0.5849625), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(allgenes, FDR<.05 & abs(logFC)>0.5849625), points(logFC, -log10(PValue), pch=20, col="green"))
with(subset(allgenes, FDR<.05 & abs(logFC)>0.5849625), textxy(logFC, -log10(PValue), cex=.6))


#pca
drug <- as.numeric(treatment)
x <- t(muscle)
pc <- prcomp(x)
plot(pc$x[,1], pc$x[,2], col=drug, main = "PCA", xlab ="PC1", ylab="PC2")
plot(pc$x[,1], pc$x[,2], col=drug, main = "PCA muscle allgenes", xlab ="PC1", ylab="PC2")


#fpkm
setwd("~/Desktop/MILES/MILES-surface/Results/")
fpkm <- read.table(file = "gene.results.fpkm.txt", header =TRUE)
dim(fpkm)
keep <- rowSums((fpkm)>1) >=2
fpkm <- fpkm[keep,]
dim(fpkm)
boxplot(log2(fpkm+1))

x <- t(fpkm)
pca <- prcomp(x)
plot(pca$x[,1], pca$x[,2], col=drug)
pca <- prcomp(log2(x+1))
plot(pca$x[,1], pca$x[,2], col=drug)









