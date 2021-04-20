---
title: "Project"
output: html_document
---

  
  The report needs to be formatted like a research paper, including title, abstract, introduction, method (analysis workflow), description of the dataset (including the cell lines studied), results (including descriptive statistics, visualization, and interpretations) and conclusion. The analyses should be centered around the questions listed in the project description. You will need report descriptive statistics, differential expression results and comparison between the two software (RNA-seq), GO term analysis (RNA-seq), differential accessibility (ATAC-seq), functional annotation (ATAC-seq) and biological interpretation. You may take a look at the data analysis section for RNA-seq or ATAC-seq data in some research papers (e.g. the references provided in the project description) to get an idea about the norm.

I include a set of slides (title only) with this email to give you a rough idea about the outline. This should not be interpreted as a template, but just basic elements to include. In the past, I didnt post this and students had been very creative. You should not be confined by this. You are encouraged to use the discussion board to discuss with others.
3. The report will be due on May 7. If you need to add more analyses or fix errors in the analyses after the presentation, feel free to do so. But the presentation is a grading component and you are expected to show most of your analyses in the presentation.

4. Someone asked a question about how to turn BigBed files into a format for running Deseq2. I obtained the script from an experience bioinformatician and put it under project module as a tab called "use BigBed files for differential analysis".  
  
  
  Questions to answer:
We are interested in the following questions:
About RNA-seq data:
1. What genes are differentially expressed across each pair of cell lines? (deseq and limma)
2. What are the functions of the genes with differential expression patterns? (onthology)
3. How consistent are the results between DEseq2 and limma voom? 
4. Construct a hierarchical tree using all the RNA-seq data, and perform clustering analysis. Describe the relationship based on your results.
About ATAC-seq data: 
5. What regions have differential chromatin patterns across each pair of cell lines? What are the genes near these regions?
6. How are differential chromatin patterns related to the expression patterns of nearby genes?
7. What are the functions of the genes with differential chromatin patterns?
8. Construct a hierarchical tree using all the ATAC-seq data and use clustering analysis to explore the pattern of cell-line specified genes. Do you get the same structure as the tree from RNA-seq data?

  

  

  
  
setwd("~/GitHub/hemopoiesis-differentiation")


f <- factor(colnames(mysample), levels=unique(colnames(mysample)))
design <- model.matrix(~0+f)
rownames(design) <- levels(f)
colSums(design)
log2fc <- as.matrix(log2fc)
fit <- lmFit(log2fc, design)



count_control1=read.csv('ENCFF247FEJ.tsv')

#remove the 'width' column
countData <- as.matrix(subset(count_control1))
#define the experimental setup 
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
#define the design formula
designFormula <- "~ group"

##Annotation



## get NCBIM37 (mm9)
ensembl67 <- useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset = "mmusculus_gene_ensembl")
## current annotation (GRCm38, mm10)
ensembl76 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

bm <- getBM(attributes=c("ensembl_gene_id", "external_gene_id", "description"), filter="ensembl_gene_id", values=as.vector(bak1), mart=ensembl67)

DESeq2Features <- data.frame(ensembl_gene_id = as.character(bak1))
DESeq2Features$ensembl_gene_id <- as.character(DESeq2Features$ensembl_gene_id)

library(dplyr)

### join them together
rowData <- dplyr::left_join(DESeq2Features, bm, by = "ensembl_gene_id")
rowData <- as(rowData, "DataFrame")
### add the annotation to the DESeq2 table
mcols(rowData(rsem.de)) <- c(mcols(rowData(rsem.de)))
#save(DESeq2Table, file = "geneCounts.RData")



filePath <- "~/GitHub/hemopoiesis-differentiation/data/"
sampleNames <- c("ENCFF858JHF", "ENCFF342WUL", "ENCFF247FEJ", "ENCFF064MKY")
countData.list <- sapply(sampleNames, function(x) read.csv(file=paste0(filePath, x, ".tsv"), header=T, sep="\t"), simplify=F)



countData.df <- do.call("cbind", countData.list)
colsToKeep <- c(1,grep("expected_count", names(countData.df)))
countss <- countData.df[,colsToKeep]
bak=countData.df[,1]
names(countss) <- c("gene_id",sampleNames)
countss[,2:4] <- round(countss[,2:4])
bak=countss[,1]
countss=na.omit(countss)
countss=cbind(Dummy=c(1:69691),countss[,-1])


   

sampleMetaData <- data.frame(cell_line=c(rep(c("HSC"), 2), rep(c("Erythroblast"),2)), Type=rep(c(rep(c("Case"), 2), rep(c("Control"), 2)),1))
rsem.in=DESeqDataSetFromMatrix(round(cts), colData = sampleMetaData, design = ~ cell_line, tidy = T)
rsem.de <- DESeq(rsem.in)










dim(rsem.de [rowSums(counts(rsem.de )) > 5, ])

results1=results(rsem.de)

results_na_o=na.omit(results1)

plotMA(results_na_o)





qvalue(results_na_o$padj,0.05)





library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)




deseq2ResDF <- as.data.frame(results_na_o)



head(deseq2ResDF)



deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .001, "Significant", NA)



ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()




ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)




which(qv$significant == TRUE)


























contrast_oe <- c("cell_line","HSC","Erythroblast")


nona=na.omit(res_tableOE)

res_tableOE <- results(rsem.de, contrast=contrast_oe, alpha = 0.05)







rsem.de <- DESeq(rsem.in)

row.names(countss)<-countss$gene_id












































y= voom(countss[,-1],mm, plot=T)



fit <- lmFit(y, mm)
head(coef(fit))


eb=eBayes(fit)

top.table=topTable(eb, n = inf)

tail(top.table)












clustering



# Transform count data using the variance stablilizing transform
deseq2VST <- vst(rsem.in)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap







# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendogram for samples
install.packages("ggdendro")
library(ggdendro)

sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))





geneid <- bak

genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
head(genes)











ATACC





openRegionBigWig <- gsub("\\.bam", "_openRegions\\.bw", sortedBAM)
openRegionRPMBigWig <- gsub("\\.bam", "_openRegionsRPM\\.bw", sortedBAM)
atacFragments_Open <- granges(atacReads_Open)
export.bw(coverage(atacFragments_Open), openRegionBigWig)






library(DESeq2)
load("ATAC_Data/ATAC_RData/countsFromATAC.RData")
metaData <- data.frame(Group, row.names = colnames(myCounts))
atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount)
atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)
plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))





library(DESeq2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tracktables)


LiverMinusHindbrain <- results(atacDDS, c("Group", "Liver", "HindBrain"), format = "GRanges")
LiverMinusHindbrain <- LiverMinusHindbrain[order(LiverMinusHindbrain$pvalue)]
LiverMinusHindbrain







library(MotifDb)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
myRes <- matchPWM(CTCF[[1]], BSgenome.Hsapiens.UCSC.hg19[["chr20"]])
toCompare <- GRanges("chr20", ranges(myRes))

read1 <- first(atacReads_Open)
read2 <- second(atacReads_Open)
Firsts <- resize(granges(read1), fix = "start", 1)
First_Pos_toCut <- shift(granges(Firsts[strand(read1) == "+"]), 4)
First_Neg_toCut <- shift(granges(Firsts[strand(read1) == "-"]), -5)


Seconds <- resize(granges(read2), fix = "start", 1)
Second_Pos_toCut <- shift(granges(Seconds[strand(read2) == "+"]), 4)
Second_Neg_toCut <- shift(granges(Seconds[strand(read2) == "-"]), -5)

test_toCut <- c(First_Pos_toCut, First_Neg_toCut, Second_Pos_toCut, Second_Neg_toCut)
cutsCoverage <- coverage(test_toCut)
cutsCoverage20 <- cutsCoverage["chr20"]
CTCF_Cuts_open <- regionPlot(cutsCoverage20, testRanges = toCompare, style = "point", 
    format = "rlelist", distanceAround = 500)
plotRegion(CTCF_Cuts_open, outliers = 0.001) + ggtitle("NucFree Cuts Centred on CTCF")















