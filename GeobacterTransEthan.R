#Geobacter Transcriptomics (mostly adapted from Damien Finn)
#install libraries. Bioconductor makes this possible. This shouldn't have to be done after the first time
biocLite(c("GenomicFeatures", "Rsamtools", "GenomicAlignments", "BiocParallel"))

#load libraries
library("GenomicFeatures")
library("Rsamtools")
library("GenomicAlignments")
library("BiocParallel")

#shortcut after saving workspace image:
load("/home/sceb/Desktop/data/users/EthanHowley/geobacterTraw/PCA_alignments_SAM_files/Wspacejan18.RData")


installed.packages() #see installed packages
#choose working directory

dir <- file.path("/home/sceb/Desktop/data/users/EthanHowley/geobacterTraw/PCA_alignments_SAM_files")
dir #check directory

list.files(dir)		#read out files in directory

#read in metadata csv
metadata <- file.path(dir, "GeobacterMeta.csv")
#read the metadata into a table from the csv
SampleTable <- read.csv(metadata, row.names = 1)
#Collect the sample names from the SampleTable into a variable
filenames <- file.path(dir, paste0(SampleTable$SampleName, ".bam"))
#check that everything exists 
file.exists(filenames)
#create list of .bam files
bamfiles <- BamFileList(filenames)

seqinfo(bamfiles) #check if this matches the assembly to count against

#set up TxDb mapping file from GFF file
gff <- file.path(dir, "GCF_000007985.2_ASM798v2_genomic.gff")
txdb <- makeTxDbFromGFF(gff, format="auto", dataSource="NCBI", circ_seqs=character()) #not sure if circ_seqs is necessary

#make a table of exons by gene
(ebg <- exonsBy(txdb, by = "gene"))
View(ebg)
seqinfo(ebg) #make sure seqinfo(ebg) matches that from the .bam files

#Create summarized overlap. This is a step where I could probably change some variables. My bam files should be paired end though

se <- summarizeOverlaps(features = ebg, reads = bamfiles, mode = "Union", singleEnd = FALSE, ignore.strand = TRUE)
#Ensure some things are actually being counted
head(assay(se), 3)
str(metadata(rowRanges(se)))

(colData(se) <- DataFrame(SampleTable))
round(colSums(assay(se))) #outputs the number of aligned sequences per sample. For Geobacter it is fairly varied but the output from the sequencer was varied so I'm not terribly concerned although maybe I should check.

#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
library(DESeq2)
dds <- DESeqDataSet(se, design = ~ description) #I don't believe I have a condition variable to control for... only comparing type
#Filter dataset to remove low count features < 20, as per Lorenz et al. 2014
dds <- dds[rowSums(counts(dds)) > 20]
#rlog transformation to make data homoskedastic
rld <- rlog(dds, blind = FALSE)
#Compare before and after transformation
par(mfrow = c( 1, 2))
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized = TRUE) [ , 1:2] + 1), pch = 16, cex = 0.3, main="before log")
plot(assay(rld)[ , 1:2], pch = 16, cex = 0.3, main="log transform")


#If happy with homoskedasticity, do some preliminary visualisations of similarity
#between samples with pheatmap
sampleDists <- dist(t(assay(rld)))
library(pheatmap)
library(RColorBrewer)
#makes a matrix with the descriptions as row names
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$description)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) ) (255)
#not really sure what plotting a heatmap with no column names tells me but it looks nice
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
#Begin DGE analysis with DESeq2 NB-GLM approach. The DESeq function will do
#sample size estimates, gene dispersion estimates and fit data to NB-GLM
dds2 <- DESeq(dds)
resdeseq <- results(dds2)
summary(resdeseq)

#plot counts
plotCounts(dds2, gene="omcZ", intgroup = "description")
plplotPCA(rld, intgroup="description")

#pretty plot count function
library(ggplot2)
prettyPlot <- function(gene1,dds){
  data <- plotCounts(dds, gene=gene1, intgroup=c("description"), returnData=TRUE)
  ggplot(data, aes(x=description, y=count, color=description ) ) +
    scale_y_log10() + 
    geom_point(position=position_jitter(width=.1,height=0), size=3) + ggtitle(gene1)
}
prettyPlot("omcZ", dds2)

#Now, look at number of genes sig different based on 0.05, 0.01 and 0.001 levels
res.05 <- results(dds2, alpha = 0.05)
resAll.05table = res.05[(res.05$padj<0.05) %in% TRUE,  ]
FBFp.05up=resAll.05table[(resAll.05table$log2FoldChange>0) %in% TRUE,]
FBFp.05down=resAll.05table[(resAll.05table$log2FoldChange<0) %in% TRUE,]

#output lists of up and down reg genes in plank v FBF
write.table(x=rownames(FBFp.05down),row.names=FALSE, col.names=FALSE, file="plankFBFdown.txt", quote=FALSE )
write.table(x=rownames(FBFp.05up),row.names=FALSE, col.names=FALSE, file="plankFBFup.txt", quote=FALSE )

summary(res.05)
res.001 <- results(dds2, alpha = 0.001)
table(res.001$padj < 0.001)
summary(res.001) #for my data, there are 389 genes with DE significant at this level. seems good?
#visualize most differential genes
library(genefilter)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 1000) #the 1000 can be changed to capture more genes, this is the top 1000 in order
#make a data frame with values normalized to mean, then a heatmap
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld) [, c("description")])
pheatmap(mat, annotation_col = df) # with 1000 genes it's hard to look at this. Let's cut it down
#same thing with 100 genes
topVarGenes100 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
mat100 <- assay(rld)[topVarGenes100, ]
mat100 <- mat100 - rowMeans(mat100)
df100 <- as.data.frame(colData(rld) [, c("description")])
pheatmap(mat100, annotation_col = df100)
#hey let's do 20, why not?
topVarGenes20 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat20 <- assay(rld)[topVarGenes20, ]
mat20 <- mat20 - rowMeans(mat20)
df20 <- as.data.frame(colData(rld) [, c("description")])
pheatmap(mat20, annotation_col = df20)
#try a different type of heatmap
library(heatmaply)
heatmaply(mat100)

#let's compare only two conditions: (nope, start from top with two conditions only)
#resHL <- results(dds2, contrast=c("description", "High", "Low")) 
resHL.05 <- results(dds2, alpha = 0.05,contrast=c("description", "High", "Low") )
table(resHL.05$padj < 0.05)
summary(resHL.05)
#heatmap top 20 high v low
#resOrderedHL <- resHL[order(resHL$padj),]
#topVarGenes20HL <- head(resOrderedHL, 20)
#selectGenesHL20 <- rownames(subset(topVarGenes20HL, ))
#mat20HL <- assay(rld)[selectGenesHL20, ]
#mat20HL <- mat20HL - rowMeans(mat20HL)
#df20HL <- as.data.frame(colData(rld) [, c("description")])
#pheatmap(mat20HL, annotation_col = df20HL)

#repeat with only high and low conditions:
metadataHL<- file.path(dir, "GeobacterMetaHL.csv")
SampleTableHL <- read.csv(metadataHL, row.names = 1)
filenamesHL <- file.path(dir, paste0(SampleTableHL$Sample_name, ".bam"))
#check that everything exists 
file.exists(filenamesHL)
#create list of .bam files
bamfilesHL <- BamFileList(filenamesHL)

seHL <- summarizeOverlaps(features = ebg, reads = bamfilesHL, mode = "Union", singleEnd = FALSE, ignore.strand = TRUE)
#Ensure some things are actually being counted
head(assay(seHL), 3)
str(metadata(rowRanges(seHL)))

(colData(seHL) <- DataFrame(SampleTableHL))
round(colSums(assay(seHL)))
ddsHL <- DESeqDataSet(seHL, design = ~ description) #I don't believe I have a condition variable to control for... only comparing type
#Filter dataset to remove low count features < 20, as per Lorenz et al. 2014
ddsHL <- ddsHL[rowSums(counts(ddsHL)) > 20]
#rlog transformation to make data homoskedastic
rldHL <- rlog(ddsHL, blind = FALSE)
#Compare before and after transformation
par(mfrow = c( 1, 2))
ddsHL <- estimateSizeFactors(ddsHL)
plot(log2(counts(ddsHL, normalized = TRUE) [ , 1:2] + 1), pch = 16, cex = 0.3, main="before log")
plot(assay(rldHL)[ , 1:2], pch = 16, cex = 0.3, main="log transform")
#If happy with homoskedasticity, do some preliminary visualisations of similarity
#between samples with pheatmap
sampleDistsHL <- dist(t(assay(rldHL)))
library(pheatmap)
library(RColorBrewer)
#makes a matrix with the descriptions as row names
sampleDistMatrixHL <- as.matrix(sampleDistsHL)
rownames(sampleDistMatrixHL) <- paste(rldHL$description)
colnames(sampleDistMatrixHL) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) ) (255)
#not really sure what plotting a heatmap with no column names tells me but it looks nice
pheatmap(sampleDistMatrixHL, clustering_distance_rows = sampleDistsHL, clustering_distance_cols = sampleDistsHL, col = colors)
#Begin DGE analysis with DESeq2 NB-GLM approach. The DESeq function will do
#sample size estimates, gene dispersion estimates and fit data to NB-GLM
ddsHL <- DESeq(ddsHL)
resHL <- results(ddsHL)
summary(resHL)
plotMA(resHL, ylim=c(-5,5), alpha=0.05, main="High vs. Low fold change, p<0.05 difference in red") 
#select points and then print names
idx <- identify(resHL$baseMean, resHL$log2FoldChange)
rownames(resHL)[idx]

##subset the differentially expressed genes padj<0.05
head(resHL[(resHL$padj<0.05) %in% TRUE,  ]) #it works!
resHL.05table = resHL[(resHL$padj<0.05) %in% TRUE,  ]
HL.05up=resHL.05table[(resHL.05table$log2FoldChange>0) %in% TRUE,]
HL.05down=resHL.05table[(resHL.05table$log2FoldChange<0) %in% TRUE,]


resHL.1table = resHL[(resHL$padj<0.1) %in% TRUE,  ]
HL.1up=resHL.1table[(resHL.1table$log2FoldChange>0) %in% TRUE,]
HL.1down=resHL.1table[(resHL.1table$log2FoldChange<0) %in% TRUE,]

#print up and down regulated gene names to txt for DAVID
write.table(x=rownames(HL.1down),row.names=FALSE, col.names=FALSE, file="HLdown.txt", quote=FALSE )
write.table(x=rownames(HL.1up),row.names=FALSE, col.names=FALSE, file="HLup.txt", quote=FALSE )

#resHLp05table<-resHL[matrix(resHL$padj < 0.05)] #this is not working
#p_1genesHL = rownames(resHL$padj<0.1) #also not working
#head(resHL.05table)
#same thing with 100 genes
topVarGenes100HL <- head(order(rowVars(assay(rldHL)), decreasing = TRUE), 100)
mat100HL <- assay(rldHL)[topVarGenes100HL, ]
mat100HL <- mat100HL - rowMeans(mat100HL)
df100HL <- as.data.frame(colData(rldHL) [, c("description")])
pheatmap(mat100HL, annotation_col = df100HL)

#how about I try some clustering with RQUBIC
library(QUBIC)
matrld <- assay(rld) #turn rlog deseq2 data into a matrix
biclustrld2 <- biclust(matrld,r=2,k=3, method=BCQU(), verbose=TRUE)
library(fields)
res <- biclust::biclust(BicatYeast, method = BCQU())
# Show the first bicluster
biclust::bicluster(matrld, res, 6)
#example: create list of bicluster 6
biclust5= biclust::bicluster(matrld, res, 5)[[1]]
#export list of gene names to a txt file
write.table(x=rownames(biclust5),row.names=FALSE, col.names=FALSE, file="biclust5genes.txt", quote=FALSE )

colnames(matrld)<- c("plank1", "plank2", "plank3", "low1", "low2", "low3", "high1","high2","high3", "fbf1", "fbf2", "fbf3")
quheatmap(matrld, biclustrld2, showlabel = TRUE, number=7 )
#make a network
net2 <- qunetwork(matrld, biclustrld, number=2, group=2, method="spearman")
if (requireNamespace("qgraph", quietly = TRUE))
  qgraph::qgraph(net2[[1]], groups = net2[[2]], layout = "spring", minimum = 0.6,
                 legend.cex = 0.5, color = c("red", "blue", "gold", "gray"), edge.label = FALSE)

#Pathway analysis
library(gage)
library(pathview)
pathview(gene.data = matrld, pathway.id="00010", species="gsu", gene.idtype = "kegg") #glycolysis


#compare only specific genes
CbCH <- c("GSU0274", "GSU3259") #(cbcl, imcH)
CbCHmat <- subset(matrld, rownames(matrld) %in% CbCH)
matCbCHdiff <- CbCHmat - rowMeans(CbCHmat)
dfCbCH <- as.data.frame(colData(rld) [, c("description")])
pheatmap(matCbCHdiff, annotation_col = dfCbCH)

#try same with only high and low
HLCbCH <- assay(rldHL)[CbCH, ]
HLCbCH <- HLCbCH - rowMeans(HLCbCH)
dfHLCbCH <- as.data.frame(colData(rldHL) [, c("description")])
pheatmap(HLCbCH, annotation_col = dfHLCbCH, cluster_cols = F)

#look at only 'cytochromes'
#download list of genes with cytochrome in description. Take symbols into GScytnames
GScytMAT <- assay(rld)[GScytNam,]
GScytMAT <- GScytMAT - rowMeans(GScytMAT)
dfGScyt <- as.data.frame(colData(rld) [, c("description")])
pheatmap(GScytMAT, annotation_col = dfGScyt, cluster_cols = F ) #pretty interesting actually

#do with just high and low
GScytMAT <- assay(rldHL)[GScytNam,]
GScytMAT <- GScytMAT - rowMeans(GScytMAT)
dfGScyt <- as.data.frame(colData(rldHL) [, c("description")])
pheatmap(GScytMAT, annotation_col = dfGScyt, cluster_cols = F )

#Can I make spatial comparisons with error bars for most interesting genes?
#choose subset of genes from rld matrix
#Create average and stdev of log2 expression counts
#make box and whisker plot?
#omcM, omcZ, omcS, omcT, omcO, omcQ, omcH, omcA


#split to different conditions
planknames <- c("plank1", "plank2", "plank3") #(planktonic samples)
plankMAT <- matrld[, colnames(matrld)%in%planknames]
Highnames <- c("high1", "high2", "high3") #(high samples)
HighMAT <- matrld[, colnames(matrld)%in%Highnames]
Lownames <- c("low1", "low2", "low3") #(low samples)
LowMAT <- matrld[, colnames(matrld)%in%Lownames]
FBFnames <- c("fbf1", "fbf2", "fbf3") #(fbf samples)
FBFMAT <- matrld[, colnames(matrld)%in%FBFnames]

#mean and SD of each condition
plankMeans <- rowMeans(plankMAT)
plankSD <- apply(plankMAT,1, sd, na.rm = TRUE)
plankMAT<- cbind(plankMAT, plankMeans)
plankMAT <- cbind(plankMAT, plankSD)

HighMeans <- rowMeans(HighMAT)
HighSD <- apply(HighMAT,1, sd, na.rm = TRUE)
HighMAT<- cbind(HighMAT, HighMeans)
HighMAT <- cbind(HighMAT, HighSD)

LowMeans <- rowMeans(LowMAT)
LowSD <- apply(LowMAT,1, sd, na.rm = TRUE)
LowMAT<- cbind(LowMAT, LowMeans)
LowMAT <- cbind(LowMAT, LowSD)

FBFMeans <- rowMeans(FBFMAT)
FBFSD <- apply(FBFMAT,1, sd, na.rm = TRUE)
FBFMAT<- cbind(FBFMAT, FBFMeans)
FBFMAT <- cbind(FBFMAT, FBFSD)


MeanSDmat <- cbind(plankMAT[, colnames(plankMAT)%in%c("plankMeans", "plankSD")], HighMAT[, colnames(HighMAT)%in%c("HighMeans", "HighSD")])
MeanSDmat <- cbind(MeanSDmat, LowMAT[, colnames(LowMAT)%in%c("LowMeans", "LowSD")], FBFMAT[, colnames(FBFMAT)%in%c("FBFMeans", "FBFSD")])

#finally do some subsetting
omc <- MeanSDmat[grep("omc", rownames(MeanSDmat)),] #subset all genes with string omc in name
write.csv(omc, file="omc-correct.csv")
#all 'pil' genes
pil <- MeanSDmat[grep("pil", rownames(MeanSDmat)),] #subset all genes with string omc in name
write.csv(pil, file="pil2.csv")


#make a heatmap between planktonic and fbf
metadataFUM<- file.path(dir, "GeobacterMetaplankFBF.csv")
SampleTableFUM <- read.csv(metadataFUM, row.names = 1)
filenamesFUM <- file.path(dir, paste0(SampleTableFUM$Sample_name, ".bam"))
#check that everything exists 
file.exists(filenamesFUM)
#create list of .bam files
bamfilesFUM <- BamFileList(filenamesFUM)

seFUM <- summarizeOverlaps(features = ebg, reads = bamfilesFUM, mode = "Union", singleEnd = FALSE, ignore.strand = TRUE)
#Ensure some things are actually being counted
head(assay(seFUM), 3)
str(metadata(rowRanges(seFUM)))

(colData(seFUM) <- DataFrame(SampleTableFUM))
round(colSums(assay(seFUM)))
ddsFUM <- DESeqDataSet(seFUM, design = ~ description) #I don't believe I have a condition variable to control for... only comparing type
#Filter dataset to remove low count features < 20, as per Lorenz et al. 2014
ddsFUM <- ddsFUM[rowSums(counts(ddsFUM)) > 20]
#rlog transformation to make data homoskedastic
rldFUM <- rlog(ddsFUM, blind = FALSE)
#Compare before and after transformation
par(mfrow = c( 1, 2))
ddsFUM <- estimateSizeFactors(ddsFUM)
plot(log2(counts(ddsFUM, normalized = TRUE) [ , 1:2] + 1), pch = 16, cex = 0.3, main="before log")
plot(assay(rldFUM)[ , 1:2], pch = 16, cex = 0.3, main="log transform")
#If happy with homoskedasticity, do some preliminary visualisations of similarity
#between samples with pheatmap
sampleDistsFUM <- dist(t(assay(rldFUM)))
library(pheatmap)
library(RColorBrewer)
#makes a matrix with the descriptions as row names
sampleDistMatrixFUM <- as.matrix(sampleDistsFUM)
rownames(sampleDistMatrixFUM) <- paste(rldHL$description)
colnames(sampleDistMatrixFUM) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) ) (255)
#not really sure what plotting a heatmap with no column names tells me but it looks nice
pheatmap(sampleDistMatrixFUM, clustering_distance_rows = sampleDistsFUM, clustering_distance_cols = sampleDistsFUM, col = colors)
#Begin DGE analysis with DESeq2 NB-GLM approach. The DESeq function will do
#sample size estimates, gene dispersion estimates and fit data to NB-GLM
ddsFUM <- DESeq(ddsFUM)
resFUM <- results(ddsFUM)
summary(resFUM)
plotMA(resFUM, ylim=c(-5,5), alpha=0.05, main="planktonic v. fbf fold change, p<0.05 difference in red") 
#same thing with 100 genes
topVarGenes100FUM <- head(order(rowVars(assay(rldFUM)), decreasing = TRUE), 100)
mat100FUM <- assay(rldFUM)[topVarGenes100FUM, ]
mat100FUM <- mat100FUM - rowMeans(mat100FUM)
df100FUM <- as.data.frame(colData(rldFUM) [, c("description")])
pheatmap(mat100FUM, annotation_col = df100FUM)

#histogram of counts b/t all 12 samples
geneCountsMean <- cbind(counts(dds), rowMeans(counts(dds)))
hist(geneCountsMean, breaks= 25, freq= TRUE, col= 'grey80', border= 'white', )


#####################
##########Counts by operon using CONDOP package##########
##############
###Need 4 inputs: DOOR ".opr" file, ".gff" file, FASTA of genome, and raw count table (se, but may need to make into matrix/DF)

#opr file downloaded from DOOR
oprG <- "D:/GeoTrans/1224.opr"
#FASTA file downloaded from NCBI
fastaG <- "D:/GeoTrans/geoPCA.fasta"
#gff file from NCBI
gffG <- "D:/GeoTrans/GCF_000007985.2_ASM798v2_genomic.gff"
#count table (hmmm, turns out I need forward and reverse info, so I may need to go back)
countG <- assay(se)
colnames(countG)<- c("plank1", "plank2", "plank3", "low1", "low2", "low3", "high1","high2","high3", "fbf1", "fbf2", "fbf3")

library(CONDOP)
#run CONDOP pre-processing
cond.in <- pre.proc(gffG, oprG, fastaG, countG, log2.expr = FALSE, sw = 100, verbose = TRUE)
