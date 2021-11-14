setwd("/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results")
# load("R/ChIPQC_DiffBind.Rdata")
# save.image("R/ChIPQC_DiffBind.Rdata")


#-----A.ChIP-seq quality assessment using ChIPQC----
## Load libraries
library("ChIPQC")
library("TxDb.Hsapiens.UCSC.hg38.knownGene") #Human referencedata
library("BiocParallel")
library("DiffBind")
library("tidyverse")
library("rtracklayer")

## Load sample data
# Read a overview sample csv file
samples <- read.csv("R/Samples.csv")
View(samples) # Look at the loaded metadata
samples <- samples[-2,]
## Create ChIPQC object
register(SerialParam()) #avoid this conflict BiocParallel
chip.Obj <- ChIPQC(samples, annotation="hg38") 

## Create ChIPQC report
ChIPQCreport(chip.Obj, reportName="ChIP QC report: INO80_ER_FOXA1", reportFolder="ChIPQCreport")
# FRiP values is the proportion of reads for that sample that overlap a peak in the consensus peakset, 
# RiP (also called FRiP) values will vary depending on the protein of interest:
# A typical good quality TF (sharp/narrow peaks) with successful enrichment would exhibit a RiP around 5% or higher.
# A good quality Pol2 (mix of sharp/narrow and dispersed/broad peaks) would exhibit a RiP of 30% or higher.
# There are also known examples of good datasets with RiP < 1% (i.e. RNAPIII or a protein that binds few sites).



#---- B. DiffBind Analysis----
#1. Input data
# 1.1 Reading in Peaksets
data.Objects <- dba(sampleSheet=samples)
# How many consensus sites were identified for this dataset
data.Objects

# 1.2 Affinity binding matrix and normalize data
# Counting reads and plot the count reads summits=250 means peaks will be 500bp, extending 250bp up and downstream of the summit)
data.count <- dba.count(data.Objects, summits=250) 
data.count
# Plot correlation PCA and heatmap using read count data to cluster the replicates and each cell line  
dba.plotHeatmap(data.count)
dba.plotPCA(data.count,  attributes=c(DBA_TISSUE, DBA_CONDITION), label=c(DBA_REPLICATE))

# normalized based on sequencing depth.
data.count <- dba.normalize(data.count)

# saveRDS(data.count, file = "R/dataCount.RDS")
# data.count <- readRDS(file = "R/dataCount.RDS")


# 2. Differential binding affinity analysis
# 2.1 Establishing contrast by group our samples based on condition we interested in
# in our case we contrast E2/Veh treatment condition
INO80.peakset <- dba.contrast(data.count, 
                              categories=DBA_TREATMENT, minMembers = 2)
                               

# DESeq2 is default method that set to used if we not change the method=DBA_EDGER, method=DESEQ2 or DBA_ALL_METHOD
# Default FDR <= 0.05 so we can adjust to 0.01
INO80.peakset <- dba.analyze(INO80.peakset) #method=DBA_DESEQ2)
#Show different results from DEG (DESEQ2)
dba.show(INO80.peakset, bContrasts=T)

#Retrieving the differentially bound sites
INO80.peakset.DB <- dba.report(INO80.peakset)

sum(INO80.peakset.DB$Fold>0)
sum(INO80.peakset.DB$Fold<0)

dba.plotVenn(INO80.peakset, contrast=1, bDB=TRUE,
             bGain=TRUE, bLoss=TRUE, bAll=FALSE)
# PROFILING AND PROFILEHEATMAP
library("profileplyr")

profiles <- dba.plotProfile(INO80.peakset, 
                            labels=list(sites=c("ER_Veh","ER_E2"))) #merge=c(DBA_TISSUE, DBA_REPLICATE), defaults merge replicates
dba.plotProfile(profiles)

# OR
profiles <- dba.plotProfile(INO80.peakset, samples=INO80.peakset$mask$ER, merge=NULL)
SpecifyingList <- GRangesList(Gain=INO80.peakset.DB[INO80.peakset.DB$Fold>0,],Loss=INO80.peakset.DB[INO80.peakset.DB$Fold<0,])
profiles <- dba.plotProfile(INO80.peakset, sites=SpecifyingList,
                            samples=list(ER_Veh=INO80.peakset$mask$ER_Veh, ER_E2=INO80.peakset$mask$ER_E2),
                            labels=list(sites=c("ER_Veh","ER_E2")), score="Fold")
dba.plotProfile(profiles)


#2.2 PLOTTING IN DIFFBIND

#Visualize which data points by using Volcano plots 
dba.plotVolcano(INO80.peakset) 
sigSites <- dba.plotVolcano(INO80.peakset, fold=log2(5)) # only highlight significant sites with at least 5x Fold Change
sigSites <- dba.plotVolcano(INO80.peakset, fold=log2(10), bLabels=TRUE)
sigSites
# Visualization using MA plots
dba.plotMA(INO80.peakset) #method=DBA_DESEQ2

# PCA samples of the different conditions cluster separately
dba.plotPCA(INO80.peakset, contrast = 1,  attributes=c(DBA_TISSUE), label=c(DBA_TISSUE)) #In contrast, 1 meaning only first contrast was generated
#Boxplots can give distribution differences between the classes 
dba.plotBox(INO80.peakset)
# display the  binding affinity correlation heat map
dba.plotHeatmap(INO80.peakset,ColAttributes = DBA_CONDITION, contrast=1) #All sample
dba.plotHeatmap(INO80.peakset, contrast=1, mask=INO80.peakset$masks$ER) #Only ESC
dba.plotHeatmap(INO80.peakset, contrast=1, mask=INO80.peakset$masks$FOXA1) #Only MEF

# Heat map based on Differential Binding sites use parameter ColAttributes = DBA_CONDITION 
Colorheatmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
dba.plotHeatmap(INO80.peakset, ColAttributes=c(DBA_TISSUE, DBA_TREATMENT), colScheme= Colorheatmap, scale="row", #contrast=1
                correlations=FALSE, method=DBA_DESEQ2, mask=INO80.peakset$masks$All) #or method=DBA_EDGER



# 2.3 Advanced occupancy analysis and overlaps

# set of sample masks are automatically associated with a DBA object in the $mask field:
names(INO80.peakset$masks)

#visualization the overlap results that as a Venn Diagram
dba.plotVenn(data.Objects, INO80.peakset$masks$ER)
dba.plotVenn(data.Objects, INO80.peakset$masks$FOXA1)

# New consensus peakset for each set of samples that share the same Tissue and Condition can be added to the DBA object using dba.peakset:
INO80_consensus <- dba.peakset(data.Objects, consensus=c(DBA_CONDITION), minOverlap=0.66) #at least represent in 2 replicates
# New DBA object can be generated consisting of only the consensus peakset
INO80_consensus <- dba(INO80_consensus, mask=INO80_consensus$masks$Consensus, minOverlap=0.66)

dba.show(INO80_consensus, INO80_consensus$masks$Consensus)

dba.plotVenn(INO80_consensus, INO80_consensus$masks$Consensus)


# overall consensus peakset, that includes peaks identified in at least two replicates of at least 1 sample group
INO80_consensus_peakset <- dba.peakset(INO80_consensus, bRetrieve=TRUE)
# This consensus peakset could then be used as the basis for the binding matrix used in dba.count:
INO80_consensus_peakset <- dba.count(NIO80.Obj, peaks=INO80_consensus_peakset)




