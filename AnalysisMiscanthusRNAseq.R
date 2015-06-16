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

DS.exp2014.h <- DESeqDataSetFromMatrix(countData = cl.exp2014[,2:97], colData = id.2014, 
                                       design = ~Harvest + Treatment*Genotype)

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

###the same analysis but now correcting for Harvest

Sys.time()
system.time(
  DS.analysis.exp2014.h <- DESeq(DS.exp2014.h, betaPrior = FALSE)
)
Sys.time()

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.h <- DESeq(DS.exp2014.h, betaPrior = FALSE, test="LRT", reduced=~Harvest + Genotype)
)
Sys.time()

resLRT.h <- results(DS.analysis.exp2014.LRT.h)

mcols(resLRT.h, use.names=TRUE)

summary(resLRT.h$padj)
sum(resLRT.h$padj<1e-04, na.rm=TRUE)

resLRT104.h <- results(DS.analysis.exp2014.LRT.h, alpha=1e-04)
sum(resLRT104.h$padj<1e-04, na.rm=TRUE)

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.i.h <- DESeq(DS.exp2014.h, betaPrior = FALSE, test="LRT", reduced=~Harvest + Treatment + Genotype)
)
Sys.time()

resLRT.i.h <- results(DS.analysis.exp2014.LRT.i.h, alpha=1e-04)
mcols(resLRT.i.h, use.names=TRUE)

summary(resLRT.i.h$padj)
sum(resLRT.i.h$padj<1e-04, na.rm=TRUE)

DS.exp2014.main.h <- DESeqDataSetFromMatrix(countData = cl.exp2014[,2:97], colData = id.2014, 
                                          design = ~ Harvest + Treatment + Genotype)

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.mT.h <- DESeq(DS.exp2014.main.h, betaPrior = FALSE, test="LRT", reduced=~Harvest + Genotype)
)
Sys.time()

resLRT.mT.h <- results(DS.analysis.exp2014.LRT.mT.h, alpha=1e-04)
mcols(resLRT.mT.h)

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.mG.h <- DESeq(DS.exp2014.main.h, betaPrior = FALSE, test="LRT", reduced=~Harvest + Treatment)
)
Sys.time()

resLRT.mG.h <- results(DS.analysis.exp2014.LRT.mG.h, alpha=1e-04)
mcols(resLRT.mG.h)
sum(resLRT.mG.h$padj<1e-04, na.rm=TRUE)

gene.names.mc <- cl.exp2014[,1]

length(gene.names.mc)

table(substr(gene.names.mc[which(resLRT.h$padj < 1e-04)], start=7, stop=9))

sum(resLRT.h$padj<1e-04, na.rm=TRUE)
sum(resLRT.i.h$padj<1e-04, na.rm=TRUE)
sum(resLRT.mT.h$padj<1e-04, na.rm=TRUE)
sum(resLRT.mG.h$padj<1e-04, na.rm=TRUE)

table(substr(gene.names.mc[which(resLRT.h$padj < 1e-04)], start=7, stop=9))
table(substr(gene.names.mc[which(resLRT.i.h$padj < 1e-04)], start=7, stop=9))
table(substr(gene.names.mc[which(resLRT.mT.h$padj < 1e-04)], start=7, stop=9))
table(substr(gene.names.mc[which(resLRT.mG.h$padj < 1e-04)], start=7, stop=9))

genes.int.h <- gene.names.mc[which(resLRT.i.h$padj < 1e-04)]
genes.int.h

dd <- data.frame(genes.int.h, round(resLRT.i.h$padj[which(resLRT.i.h$padj < 1e-04)], digits=20))

dd.order <- dd[order(dd[,2]),]

xtable(dd.order, digits=10)

w.i <- which(resLRT.i.h$padj < 1e-04)
w.m <- which(resLRT.mT.h$padj < 1e-04) 
w.i%in%w.m
w.i[w.i%in%w.m]
w.m[w.m%in%w.i]

genes.mT.h <- gene.names.mc[w.m] 
indi <- which(genes.mT.h%in%genes.int.h)
genes.mT.h.pure <- genes.mT.h[-indi]

write.table(genes.mT.h.pure, "DEgenes_mainEffect_Treatment_alphabetically.txt", sep="\t")

write.table(dd.order, "DEgenes_interaction_byP-Value.txt", sep="\t")

#just for fun: check for differential expression for harvest

Sys.time()
system.time(
  DS.analysis.exp2014.LRT.h.h <- DESeq(DS.exp2014.h, betaPrior = FALSE, test="LRT", reduced=~Treatment*Genotype)
)
Sys.time()
resLRT.h.h <- results(DS.analysis.exp2014.LRT.h.h, alpha=1e-04)
mcols(resLRT.h.h)
sum(resLRT.h.h$padj<1e-04, na.rm=TRUE)

sessionInfo()



# > sessionInfo()
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
#   [1] MASS_7.3-40               DESeq2_1.8.1              RcppArmadillo_0.5.100.1.0 Rcpp_0.11.6               GenomicRanges_1.20.3     
# [6] GenomeInfoDb_1.4.0        IRanges_2.2.1             S4Vectors_0.6.0           BiocGenerics_0.14.0       rgl_0.95.1247            
# [11] gridExtra_0.9.1           ggplot2_1.0.1             gplots_2.17.0             RColorBrewer_1.1-2        lattice_0.20-31          
# [16] XLConnect_0.2-11          XLConnectJars_0.2-9      
# 
# loaded via a namespace (and not attached):
#   [1] genefilter_1.50.0    gtools_3.4.2         locfit_1.5-9.1       reshape2_1.4.1       splines_3.2.0        rJava_0.9-6         
# [7] colorspace_1.2-6     survival_2.38-1      XML_3.98-1.1         foreign_0.8-63       DBI_0.3.1            BiocParallel_1.2.1  
# [13] lambda.r_1.1.7       plyr_1.8.2           stringr_1.0.0        munsell_0.4.2        gtable_0.1.2         futile.logger_1.4.1 
# [19] caTools_1.17.1       labeling_0.3         latticeExtra_0.6-26  Biobase_2.28.0       geneplotter_1.46.0   AnnotationDbi_1.30.1
# [25] proto_0.3-10         acepack_1.3-3.3      KernSmooth_2.23-14   xtable_1.7-4         scales_0.2.4         gdata_2.16.1        
# [31] Hmisc_3.16-0         annotate_1.46.0      XVector_0.8.0        digest_0.6.8         stringi_0.4-1        tools_3.2.0         
# [37] bitops_1.0-6         magrittr_1.5         RSQLite_1.0.0        Formula_1.2-1        cluster_2.0.1        futile.options_1.0.0
# [43] rpart_4.1-9          nnet_7.3-9          

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

pv.harvest <- rep(NA, dim(cl.exp2014)[1])

system.time(
for(i in 1:dim(cl.exp2014)[1]){
  y <- cl.exp2014[i,2:97]
  mod <- glm.nb(unlist(y) ~ Treatment*Genotype + Harvest, data=id.2014) 
  pv.harvest[i] <- summary(aov(mod))[[1]]$Pr[3] #p-value from the anova table for Harvest
  print(i)  
}
)

pv.block <- rep(NA, dim(cl.exp2014)[1])
pv.row <- rep(NA, dim(cl.exp2014)[1])
pv.col <- rep(NA, dim(cl.exp2014)[1])


system.time(
  for(i in 1:dim(cl.exp2014)[1]){
    y <- cl.exp2014[i,2:97]
    mod1 <- glm.nb(unlist(y) ~ Treatment*Genotype + block, data=id.2014) 
    pv.block[i] <- summary(aov(mod1))[[1]]$Pr[3] #p-value from the anova table for Harvest
    print(i)  
  }
)

summary(pv.harvest)
length(which(pv.harvest<0.05))
length(which(pv.harvest<0.01))
length(which(pv.harvest<0.001))
length(which(pv.harvest<0.0001))

summary(pv.block)
length(which(pv.block<0.05))
length(which(pv.block<0.01))
length(which(pv.block<0.001))
length(which(pv.block<0.0001))

which(resLRT.i$padj < 1e-04)%in%which(pv.harvest<0.0001)
which(resLRT.i$padj < 1e-04)%in%which(pv.harvest<0.001)
which(resLRT.i$padj < 1e-04)%in%which(pv.harvest<0.01)
which(resLRT.i$padj < 1e-04)%in%which(pv.harvest<0.05)

sum(which(resLRT$padj < 1e-04)%in%which(pv.harvest<0.0001))
sum(which(resLRT.mT$padj < 1e-04)%in%which(pv.harvest<0.0001))

which(resLRT.i$padj < 1e-04)%in%which(pv.harvest<0.001)
which(resLRT.i$padj < 1e-04)%in%which(pv.harvest<0.01)
which(resLRT.i$padj < 1e-04)%in%which(pv.harvest<0.05)


#error message with row/column
# 
# system.time(
#   for(i in 1:dim(cl.exp2014)[1]){
#     y <- cl.exp2014[i,2:97]
#     mod2 <- glm.nb(unlist(y) ~ Treatment*Genotype + newrow14 + newcol14, data=id.2014) 
#     pv.row[i] <- summary(aov(mod2))[[1]]$Pr[3] #p-value from the anova table for Harvest
#     pv.col[i] <- summary(aov(mod2))[[1]]$Pr[4] #p-value from the anova table for Harvest
#     print(i)  
#   }
# )









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





