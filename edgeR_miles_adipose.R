######AMEYA KULKARNI#############
######MILES RNASEQ ANALYSIS######
######ADIPOSE########

setwd("Desktop/MILES/")

##Read mus-matrix file created from rsem
fat <- read.table(file="fat-matrix.txt", header=TRUE)

fatnames <- colnames(fat)
fatnames <- strsplit(fatnames, split = "[.]")
fatnames <- lapply(fatnames, function(x){paste(x[1],x[2],x[3],x[4], sep=".")})
fatnames<-fatnames[-30]
unlist(fatnames)


sfat <- strsplit(as.character(fatnames), split="[.]")
sfat <- lapply(sfat, function(x){x[2]})
unlist(sfat)

tfat <- strsplit(as.character(fatnames), split="[.]")
tfat <- lapply(tfat, function(x){x[4]})
unlist(tfat)

vfat <- strsplit(as.character(fatnames), split="[.]")
vfat <- lapply(vfat, function(x){x[3]})
unlist(vfat)

tarfat <- cbind(fatnames, sfat, tfat, vfat)
tarfat <- data.frame(tarfat)

#Generating groups for edgeR
tfat <- factor(as.character(tarfat$tfat))
sfat <- factor(as.character(tarfat$sfat))
vfat <- factor(as.character(tarfat$vfat))

library(edgeR)

#Round off expected counts to the nearest lower integer
fat <- round(fat, digits=0)

#Data matrix for expected counts without columns containing subjects with one biopsy
fat <- fat[-30] #Use this counts matrix for all further analysis

#Generating a DGEList specific for edgeR analysis in R
genefat <- DGEList(counts=fat, group=tfat)
genefat <- calcNormFactors(genefat)
cpmfat <- cpm(genefat$counts)
logcpmfat <- log(cpm(genefat$counts))

#Filtering
keepfat <- rowSums((cpm(genefat))>2) >=2
genefat <- genefat[keepfat, , keep.lib.sizes=FALSE]

#Design restructuring

designfat <- model.matrix(~tfat+sfat+vfat)
yfat <- estimateGLMCommonDisp(genefat,designfat)
yfat <- estimateGLMTrendedDisp(yfat,designfat)
yfat <- estimateGLMTagwiseDisp(yfat,designfat)


fitfat <- glmFit(yfat,designfat)
lrtfat <- glmLRT(fitfat, coef=2)
topTags(lrtfat)

plotBCV(yfat)
ofat <- order(lrtfat$table$PValue)
cpm(yfat)[o[1:10],]
summary(defat <- decideTestsDGE(lrtfat))
detagsfat <- rownames(yfat)[as.logical(defat)]
plotSmear(lrtfat, de.tags=detagsfat)
abline(h=c(-2, 2), col="green")

genenamesfat <- strsplit(as.character(rownames((topTags(lrtfat, n= Inf)))), split="[_]")
genenamesfat <- lapply(genenamesfat, function(x){x[2]})

fatdegtsv <- topTags(lrtfat, n=Inf)
genenamesfat <- data.frame(unlist(genenamesfat))
fatdegtsv <- cbind(genenamesfat, fatdegtsv)
write.csv (fatdegtsv, file="miles_fat_deg.csv")




#Volcano plots
allgenesfat <- read.delim(file="allgenesfat.txt")

#For FC>2 or <0.5
with(allgenesfat, plot(logFC, -log10(PValue), pch=20, main="DEG in adipose after Metformin treatment", xlim=c(-5.5,6)))
with(subset(allgenesfat, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(allgenesfat, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(allgenesfat, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
with(subset(allgenesfat, FDR<.05 & abs(logFC)>1), textxy(logFC, -log10(PValue), labs=fatgenesym, cex=.6))

#for FC>1.5 or FC <(1/1.5)
with(allgenesfat, plot(logFC, -log10(PValue), pch=20, main="DEG in muscle after Metformin treatment", xlim=c(-4,5)))
with(subset(allgenesfat, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(allgenesfat, abs(logFC)>0.5849625), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(allgenesfat, FDR<.05 & abs(logFC)>0.5849625), points(logFC, -log10(PValue), pch=20, col="green"))
with(subset(allgenesfat, FDR<.05 & abs(logFC)>0.5849625), textxy(logFC, -log10(PValue), cex=.6))

fatgenesym <- unlist(lapply(strsplit(as.character(allgenesfat$GENCODE.ID), split = "[_]"), function(x){x[2]}))
allgenesfat <- cbind(fatgenesym, allgenesfat)