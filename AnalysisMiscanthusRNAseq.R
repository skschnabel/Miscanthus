#Analysis for the 2014 data only

#read in the miscanthus data and store into a workspace
setwd("D:/Ausgelagert/WATBIO/Miscanthus_Counts/")

exp2014 <- read.table("D:/Ausgelagert/WATBIO/Miscanthus_Counts/htseq_T1.2-Miscanthus-drought-Exp2014_batch2_bigTable.tsv/htseq_T1.2-Miscanthus-drought-Exp2014_batch2_bigTable.tsv", 
                      sep="\t", header=TRUE)
library(XLConnect)
library(lattice)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(gridExtra)
library(rgl)
library(DESeq2)


#meta-data
meta.wb <- loadWorkbook("Additional Information_Marta_new.xlsx")
meta.data <- readWorksheet(meta.wb, sheet = "info", header = TRUE)

meta.data.2014 <- meta.data[meta.data$Experiment==2014,]
meta.data.2014 <- meta.data.2014[order(meta.data.2014$numnum), ]

counts.sq.2014 <- prcomp(t(sqrt(exp2014[,2:97])), center=TRUE, .scale=TRUE)
qplot(counts.sq.2014$x[,1], counts.sq.2014$x[,2], col=meta.data.2014$Genotype, size=I(5), pch=meta.data.2014$Treatment)

plot3d(counts.sq.2014$x[,1], counts.sq.2014$x[,2], counts.sq.2014$x[,3], col=as.numeric(as.factor(meta.data.2014$Genotype)), size=15, 
       pch=as.numeric(meta.data.2014$Treatment)+15) 

#####################

####PLOTS FIELD DESIGN
##############
###thinking about the field design

colnames(meta.data.2014)
ii <- as.numeric(as.character(interaction(meta.data.2014$Genotype, meta.data.2014$Treatment)))

mm <- matrix(ncol=25, nrow=7, NA)

meta.data.2014$newrow14 <- as.numeric(as.character(meta.data.2014$newrow14))
meta.data.2014$newcol14 <- as.numeric(as.character(meta.data.2014$newcol14))

mm[cbind(meta.data.2014$newrow14, meta.data.2014$newcol14)] <- ii

image(t(mm[nrow(mm):1,]), x=1:25, y=-7:-1)

ifi <- as.factor(ii)
library(RColorBrewer)
coco <- brewer.pal(10, "Paired")

mm[cbind(meta.data.2014$newrow14, meta.data.2014$newcol14)] <- ifi

image(t(mm[nrow(mm):1,]), x=1:25, y=-7:-1, col=coco)

ccc <- cbind((-1)*meta.data.2014$newrow14, meta.data.2014$newcol14)
text(x=ccc[,2], y=ccc[,1], labels=ii)


################
###clean data
###see also functions at the end of code

m10.2014 <- which(apply(exp2014[,2:97], 1, mean.gr.0)>10)
cl.exp2014 <- exp2014[m10.2014,]

##################

head(cl.exp2014)

cl.sqrt.counts.2014 <- sqrt(cl.exp2014[,2:97])
head(cl.sqrt.counts.2014)

id.2014 <-meta.data.2014
id.2014$Genotype <- as.factor(id.2014$Genotype)
id.2014$Treatment <- as.factor(id.2014$Treatment)
id.2014$Harvest <- as.factor(id.2014$Harvest)
id.2014$newrow14 <- as.factor(id.2014$newrow14)
id.2014$newcol14 <- as.factor(id.2014$newcol14)

DS.exp2014 <- DESeqDataSetFromMatrix(countData = cl.exp2014[,2:97], colData = id.2014, 
                                     design = ~Treatment*Genotype)

Sys.time()
system.time(
  DS.analysis.exp2014 <- DESeq(DS.exp2014, betaPrior = FALSE)
)
Sys.time()

Sys.time()
system.time(
  DS.analysis.exp2014.LRT <- DESeq(DS.exp2014, betaPrior = FALSE, test="LRT", reduced=~Genotype)
)
Sys.time()

resLRT <- results(DS.analysis.exp2014.LRT)

mcols(resLRT, use.names=TRUE)

summary(resLRT$padj)
sum(resLRT$padj<1e-04, na.rm=TRUE)

resLRT104 <- results(DS.analysis.exp2014.LRT, alpha=1e-04)
sum(resLRT104$padj<1e-04, na.rm=TRUE)

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.i <- DESeq(DS.exp2014, betaPrior = FALSE, test="LRT", reduced=~Treatment + Genotype)
)
Sys.time()

resLRT.i <- results(DS.analysis.exp2014.LRT.i, alpha=1e-04)


DS.exp2014.main <- DESeqDataSetFromMatrix(countData = cl.exp2014[,2:97], colData = id.2014, 
                                           design = ~ Treatment + Genotype)

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.mT <- DESeq(DS.exp2014.main, betaPrior = FALSE, test="LRT", reduced=~Genotype)
)
Sys.time()

resLRT.mT <- results(DS.analysis.exp2014.LRT.mT, alpha=1e-04)
mcols(resLRT.mT)

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.mG <- DESeq(DS.exp2014.main, betaPrior = FALSE, test="LRT", reduced=~Treatment)
)
Sys.time()

resLRT.mG <- results(DS.analysis.exp2014.LRT.mG, alpha=1e-04)
mcols(resLRT.mG)
sum(resLRT.mG$padj<1e-04, na.rm=TRUE)

gene.names.mc <- cl.exp2014[,1]

length(gene.names.mc)

table(substr(gene.names.mc[which(resLRT$padj < 1e-04)], start=7, stop=9))

sum(resLRT$padj<1e-04, na.rm=TRUE)
sum(resLRT.i$padj<1e-04, na.rm=TRUE)
sum(resLRT.mT$padj<1e-04, na.rm=TRUE)
sum(resLRT.mG$padj<1e-04, na.rm=TRUE)

table(substr(gene.names.mc[which(resLRT$padj < 1e-04)], start=7, stop=9))
table(substr(gene.names.mc[which(resLRT.i$padj < 1e-04)], start=7, stop=9))
table(substr(gene.names.mc[which(resLRT.mT$padj < 1e-04)], start=7, stop=9))
table(substr(gene.names.mc[which(resLRT.mG$padj < 1e-04)], start=7, stop=9))

genes.int <- gene.names.mc[which(resLRT.i$padj < 1e-04)]
genes.int



# sl <- slot(resLRT.i, "listData")
# sl.m <- matrix(unlist(sl), nrow=16343, ncol=6)
# dim(sl.m)
# 
# sl.m[which(resLRT.i$padj < 1e-04),]
# 
# levelplot(cor(cl.exp2014[,2:97]), xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
# 
# which(id.2014$Harv==1)
# 
# seqqi <- 2:97
# seqqi1 <- seqqi[which(id.2014$Harv==1)]
# seqqi2 <- seqqi[which(id.2014$Harv==2)]
# 
# levelplot(cor(cl.exp2014[,seqqi1]), xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
# levelplot(cor(cl.exp2014[,seqqi2]), xlab=NULL, ylab=NULL, scales=list(draw=FALSE))


y <- cl.exp2014[1,2:97]

library(MASS)

m1 <- glm.nb(unlist(y) ~ Treatment*Genotype + Harvest, data=id.2014) 
summary(m1)
summary(aov(m1))
summary(aov(m1))[[1]]$Pr[3] #p-value from the anova table for Harvest










##################
min.gr.0 <- function(x){
  ifelse(sum(x)>0, min(x[x>0]), -999)
}

mean.gr.0 <- function(x){
  ifelse(sum(x)>0, mean(x[x>0]), -9999)
}

number.gr.0 <- function(x){
  ifelse(sum(x)>0, length(x[x>0]), 0)
}

number.mean.gr.10 <- function(x){
  ifelse(sum(x)>0, length(x[mean(x[x>0])]), 0)
}





