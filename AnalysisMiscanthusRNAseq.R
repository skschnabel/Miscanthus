#first analyses of the miscanthus data -- meta data and RNAseq

# R version 3.2.0 (2015-04-16)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
# [4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252    
# 
# attached base packages:
#   [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DESeq2_1.8.1              RcppArmadillo_0.5.100.1.0 Rcpp_0.11.6               GenomicRanges_1.20.3      GenomeInfoDb_1.4.0       
# [6] IRanges_2.2.1             S4Vectors_0.6.0           BiocGenerics_0.14.0       rgl_0.95.1247             gridExtra_0.9.1          
# [11] ggplot2_1.0.1             gplots_2.17.0             RColorBrewer_1.1-2        lattice_0.20-31           XLConnect_0.2-11         
# [16] XLConnectJars_0.2-9      
# 
# loaded via a namespace (and not attached):
#   [1] genefilter_1.50.0    gtools_3.4.2         locfit_1.5-9.1       reshape2_1.4.1       splines_3.2.0        rJava_0.9-6         
# [7] colorspace_1.2-6     survival_2.38-1      XML_3.98-1.1         foreign_0.8-63       DBI_0.3.1            BiocParallel_1.2.1  
# [13] lambda.r_1.1.7       plyr_1.8.2           stringr_1.0.0        munsell_0.4.2        gtable_0.1.2         futile.logger_1.4.1 
# [19] caTools_1.17.1       labeling_0.3         latticeExtra_0.6-26  Biobase_2.28.0       geneplotter_1.46.0   AnnotationDbi_1.30.1
# [25] proto_0.3-10         acepack_1.3-3.3      KernSmooth_2.23-14   xtable_1.7-4         scales_0.2.4         gdata_2.16.1        
# [31] Hmisc_3.16-0         annotate_1.46.0      XVector_0.8.0        digest_0.6.8         stringi_0.4-1        tools_3.2.0         
# [37] bitops_1.0-6         magrittr_1.5         RSQLite_1.0.0        Formula_1.2-1        cluster_2.0.1        futile.options_1.0.0
# [43] MASS_7.3-40          rpart_4.1-9          nnet_7.3-9    

#read in the miscanthus data and store into a workspace
#updated file 2013
setwd("D:/Ausgelagert/WATBIO/Miscanthus_Counts/")

exp2013 <- read.table("D:/Ausgelagert/WATBIO/Miscanthus_Counts/htseq_T1.2-Miscanthus-drought-heat-Exp2013_batch3_bigTable_new.tsv/htseq_T1.2-Miscanthus-drought-heat-Exp2013_batch3_bigTable.tsv",
                      sep="\t", header=TRUE)

exp2014 <- read.table("D:/Ausgelagert/WATBIO/Miscanthus_Counts/htseq_T1.2-Miscanthus-drought-Exp2014_batch2_bigTable.tsv/htseq_T1.2-Miscanthus-drought-Exp2014_batch2_bigTable.tsv", 
                      sep="\t", header=TRUE)


#meta-data
library(XLConnect)

meta.wb <- loadWorkbook("Additional Information_Marta_new.xlsx")

meta.data <- readWorksheet(meta.wb, sheet = "info", header = TRUE)

meta.data.2013 <- meta.data[meta.data$Experiment==2013,]
meta.data.2013 <- meta.data.2013[order(meta.data.2013$numnum), ]

meta.data.2014 <- meta.data[meta.data$Experiment==2014,]
meta.data.2014 <- meta.data.2014[order(meta.data.2014$numnum), ]


library(lattice)
library(RColorBrewer)
library(gplots)
library(ggplot2)

counts.sq.2013 <- prcomp(t(sqrt(exp2013[,2:91])), center=TRUE, .scale=TRUE)
p2013 <- qplot(counts.sq.2013$x[,1], counts.sq.2013$x[,2], col=meta.data.2013$Genotype, size=I(5), pch=meta.data.2013$Treatment)

counts.sq.2014 <- prcomp(t(sqrt(exp2014[,2:97])), center=TRUE, .scale=TRUE)
p2014 <- qplot(counts.sq.2014$x[,1], counts.sq.2014$x[,2], col=meta.data.2014$Genotype, size=I(5), pch=meta.data.2014$Treatment)

library(gridExtra)

grid.arrange(p2013, p2014, ncol=2)

#pdf("PCA_all_sqrt_counts_24032015.pdf", width=12, height=4)
grid.arrange(p2013, p2014, ncol=2)
#dev.off()

library(rgl)
plot3d(counts.sq.2014$x[,1], counts.sq.2014$x[,2], counts.sq.2014$x[,3], col=as.numeric(as.factor(meta.data.2014$Genotype)), size=15, 
       pch=as.numeric(meta.data.2014$Treatment)+15) 


################
###clean data

summary(apply(exp2013[,2:91], 2, min))
summary(apply(exp2013[,2:91], 2, max))
summary(apply(exp2013[,2:91], 2, sum))

summary(apply(exp2013[,2:91], 1, mean))

table(apply(exp2013[,2:91], 1, sum))[1:10]
#we have 33032 genes in total
#5730 genes have no counts at all

table(apply(exp2013[,2:91], 1, mean.gr.0))[1:20]

#keeping in mind that we have about 90 samples

length(which((apply(exp2013[,2:91], 1, mean.gr.0)>10)))
#16674 genes have a mean count larger than 10 (only based on the non-zero counts)
length(which((apply(exp2013[,2:91], 1, mean)>10)))
#16438 genes have a mean count larger than 10 


length(which((apply(exp2013[,2:91], 1, sum)>10)))
#24233 genes have a library size larger than 10

length(which((apply(exp2013[,2:91], 1, sum)>50)))
#21885 genes have a library size larger than 50

m10.2013 <- which(apply(exp2013[,2:91], 1, mean.gr.0)>10)

cl.exp2013 <- exp2013[m10.2013,]

m10.2014 <- which(apply(exp2014[,2:97], 1, mean.gr.0)>10)

cl.exp2014 <- exp2014[m10.2014,]

##################

library(DESeq2)

head(cl.exp2014)

cl.sqrt.counts.2014 <- sqrt(cl.exp2014[,2:97])
head(cl.sqrt.counts.2014)

id.2014 <-meta.data.2014
id.2014$Genotype <- as.factor(id.2014$Genotype)
id.2014$Treatment <- as.factor(id.2014$Treatment)
id.2014$Harvest <- as.factor(id.2014$Harvest)
id.2014$row <- as.factor(id.2014$row)
id.2014$column <- as.factor(id.2014$column)


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

# 
# mcols(r1)
# 
# r2 <- results(DS.analysis.exp2014, name="Treamtent_2_vs_1")
# r2 <- results(DS.analysis.exp2014, name="Treatment_2_vs_1")
# mcols(r2)
# summary(r2$padj)
# sum(r2$padj<0.01, na.rm=TRUE)
# replications(meta.data.2014)
# replications(meta.data.2013)
# table(meta.data.2013)
# table(meta.data.2013$Treat)
# table(meta.data.2013$Harvest)

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

as.matrix(resLRT.i@listData)[which(resLRT.i$padj < 1e-04),]

sl <- slot(resLRT.i, "listData")
sl.m <- matrix(unlist(sl), nrow=16343, ncol=6)
dim(sl.m)

sl.m[which(resLRT.i$padj < 1e-04),]


levelplot(cor(cl.exp2014[,2:97]), xlab=NULL, ylab=NULL, scales=list(draw=FALSE))

which(id.2014$Harv==1)

seqqi <- 2:97
seqqi1 <- seqqi[which(id.2014$Harv==1)]
seqqi2 <- seqqi[which(id.2014$Harv==2)]

levelplot(cor(cl.exp2014[,seqqi1]), xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
levelplot(cor(cl.exp2014[,seqqi2]), xlab=NULL, ylab=NULL, scales=list(draw=FALSE))


y <- cl.exp2014[1,2:97]

library(MASS)

m1 <- glm.nb(unlist(y) ~ Treatment*Genotype + Harvest, data=id.2014) 
summary(m1)
summary(aov(m1))
summary(aov(m1))[[1]]$Pr[3] #p-value from the anova table for Harvest






#####################

####PLOTS FIELD DESIGN
##############
###thinking about the field design

colnames(meta.data.2014)
ii <- as.numeric(as.character(interaction(meta.data.2014$Genotype, meta.data.2014$Treatment)))

mm <- matrix(ncol=17, nrow=6, NA)

meta.data.2014$row <- as.numeric(as.character(meta.data.2014$row))
meta.data.2014$column <- as.numeric(as.character(meta.data.2014$column))

mm[cbind(meta.data.2014$row, meta.data.2014$column)] <- ii

image(t(mm[nrow(mm):1,]), x=1:17, y=-6:-1)

ifi <- as.factor(ii)
library(RColorBrewer)
coco <- brewer.pal(10, "Paired")

mm[cbind(meta.data.2014$row, meta.data.2014$column)] <- ifi

image(t(mm[nrow(mm):1,]), x=1:17, y=-6:-1, col=coco)


ccc <- cbind((-1)*meta.data.2014$row, meta.data.2014$column)

text(x=ccc[,2], y=ccc[,1], labels=ii)



#################

iii <- as.numeric(as.character(interaction(meta.data.2013$Genotype, meta.data.2013$Treatment)))

mmm <- matrix(ncol=14, nrow=10, NA)

meta.data.2013$row <- as.numeric(as.character(meta.data.2013$row))
meta.data.2013$column <- as.numeric(as.character(meta.data.2013$column))

mmm[cbind(meta.data.2013$row, meta.data.2013$column)] <- iii

image(t(mmm[nrow(mmm):1,]), x=1:14, y=-10:-1)

iffi <- as.factor(iii)
mmm[cbind(meta.data.2013$row, meta.data.2013$column)] <- iffi

image(t(mmm[nrow(mmm):1,]), x=1:14, y=-10:-1, col=coco)

cccc <- cbind((-1)*meta.data.2013$row, meta.data.2013$column)

text(x=cccc[,2], y=cccc[,1], labels=iii)

pdf("Designs-Miscanthus.pdf", width=10, height=5)
par(mfrow=c(1,2))

image(t(mmm[nrow(mmm):1,]), x=1:14, y=-10:-1, col=coco, main="Design 2013", xlab="", ylab="")
cccc <- cbind((-1)*meta.data.2013$row, meta.data.2013$column)
text(x=cccc[,2], y=cccc[,1], labels=iii, cex=0.75)

image(t(mm[nrow(mm):1,]), x=1:17, y=-6:-1, col=coco, main="Design 2014", xlab="", ylab="")
ccc <- cbind((-1)*meta.data.2014$row, meta.data.2014$column)
text(x=ccc[,2], y=ccc[,1], labels=ii, cex=0.75)
dev.off()


#################







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





# summary.soton[,1] <- apply(soton.gxd.data, 1, min)
# summary.soton[,2] <- apply(soton.gxd.data, 1, min.gr.0)
# summary.soton[,3] <- apply(soton.gxd.data, 1, max)
# summary.soton[,4] <- apply(soton.gxd.data, 1, sum)
# summary.soton[,5] <- apply(soton.gxd.data, 1, mean)
# summary.soton[,6] <- apply(soton.gxd.data, 1, mean.gr.0)
# summary.soton[,7] <- apply(soton.gxd.data, 1, number.gr.0)
