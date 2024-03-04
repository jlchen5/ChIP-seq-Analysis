#下载安装ChIPseeker注释相关的包
source ("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")

#下载人的基因和lincRNA的TxDb对象
biocLite("org.Dm.eg.db")
biocLite("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
biocLite("clusterProfiler")

#载入各种包
library("ChIPseeker")
library("clusterProfiler")
library("org.Dm.eg.db")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

#读入peak文件
orc2 <- readPeakFile("./idr_out.bed/nanog_idr-bed")
zld <- readPeakFile("./idr_out.bed/pou5f1_idr-bed")

# 制作多个样本比较的list
peaks <- list(Nanog=nanog,Pou5f1=pou5f1)
# promotor区间范围可以自己设定
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
#annotatePeak传入annoDb参数,可进行基因ID转换（Entrez，ENSEMBL，SYMBOL，GENENAME）
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000),verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")

#可视化
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList,title="Distribution of transcription factor-binding loci \n relative to TSS")