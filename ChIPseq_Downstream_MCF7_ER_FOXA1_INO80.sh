bowtie2Ref=/media/hkh/8TB/XUANTHANG/References/Refrence_Human_Genome/Homo_sapiens_UCSC_hg38/Sequence/Bowtie2Index/genome
Reference=/media/hkh/8TB/XUANTHANG/References/Refrence_Human_Genome
Align="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/sambamba"
MACS2="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/MACS2"
log="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/MACS2/MACS2log"
bigWig="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/bigWig"
bedGraphToBigWig="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/bedGraphToBigWig"
intervene="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/intervene"
bedtools="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/bedtools"

# DOWNSTREAM ANALYSIS 
################################################################################################################################################
# 5 deepTools Visualization

# 5.1 Generate bigWig file from BAM file
# Using bedGraphToBigWig tool (with treat.pileup.bdg) or bamCompare/bamCoverage from bedTools
mkdir -p $Results/bigWig/bigWiglogs $Results/bedGraphToBigWig

#Create the Chrome size what using for alignment (fasta file)
samtools faidx $bowtie2Ref.fa | cut -f 1,2 $bowtie2Ref.fa.fai > $MACS2/chrom.sizes

# Create bigWig file from bedGraph file using bedGraphToBigWig file tool
#for i in $MACS2/*_treat_pileup.bdg ; do
#	echo $i
#	base=`basename -s _treat_pileup.bdg $i` 
#	sort -k1,1 -k2,2n $i > $MACS2/${base}_treat_pileup_sorted.bdg
#	bedGraphToBigWig $MACS2/${base}_treat_pileup_sorted.bdg $MACS2/chrom.sizes $bedGraphToBigWig/${base}.bw
#done

#OR
### Convert .bam files to .bigWig(.bw) files
for i in $Align/*.bam ; do
	echo $i;
	base=`basename -s .bam $i`
	bamCoverage -p 6 \
                --centerReads --extendReads 150  --binSize 10 \
                --normalizeUsing RPGC --effectiveGenomeSize 3209286105 \
                --ignoreForNormalization chrX chrY \
                -b $Align/${base}.bam \
                -o $bigWig/${base}.bw 2> $bigWig/bigWiglogs/${base}_bamCoverage.log
done

#Merge bigWig replicates
# Merge three files for ER, FOXA1, and INO80.
# conda install -c bioconda ucsc-bigwigmerge
bigWigMerge $bigWig/SRR1635435.bw $bigWig/SRR1635436.bw $bigWig/Input_Veh.bedgraph
bigWigMerge $bigWig/SRR1635437.bw $bigWig/SRR1635438.bw $bigWig/Input_E2.bedgraph

bigWigMerge $bigWig/SRR1635443.bw $bigWig/SRR1635444.bw $bigWig/ER_Veh.bedgraph
bigWigMerge $bigWig/SRR1635445.bw $bigWig/SRR1635446.bw $bigWig/ER_E2.bedgraph
bigWigMerge $bigWig/SRR1635459.bw $bigWig/SRR1635460.bw $bigWig/FOXA1_Veh.bedgraph
bigWigMerge $bigWig/SRR1635461.bw $bigWig/SRR1635462.bw $bigWig/FOXA1_E2.bedgraph


# Sort and get average scores.
for i in $bigWig/*.bedgraph ; do
    echo $i;
    base=`basename -s .bedgraph $i`
    sort -k1,1 -k2,2n $bigWig/${base}.bedgraph | \
    awk 'BEGIN{OFS = "\t"}($5 = $4/3){print $1,$2,$3,$5}' > $bigWig/${base}_adjusted.bedgraph 
    bedGraphToBigWig $bigWig/${base}_adjusted.bedgraph $MACS2/chrom.sizes $bedGraphToBigWig/${base}.bw
done







#5.2 Working with MACS2 Peakcalling results  
# Obtain intersection of the replicates for each sample.
mkdir -p $Results/bedtools $Results/intervene

for i in $MACS2/*_rep1_peaks.narrowPeak ; do
    echo $i;
    base=`basename -s _rep1_peaks.narrowPeak $i`
    intervene venn -i $MACS2/${base}_rep1_peaks.narrowPeak $MACS2/${base}_rep2_peaks.narrowPeak \
                   --save-overlaps -o $intervene/${base} \
                   --names=${base}_rep1,${base}_rep2
done

## Concatenate the peaks that are in overlapped regions.

#If using more than 2 replicates,choose Peaks that are common to at least two out of three replicates are selected.
#cat 011_rep2_rep3.bed 101_rep1_rep3.bed 110_rep1_rep2.bed 111_rep1_rep2_rep3.bed > WT_rep.bed 

for i in $MACS2/*_rep1_peaks.narrowPeak ; do
    echo $i;
    base=`basename -s _rep1_peaks.narrowPeak $i`
    cat $intervene/${base}/sets/11_*.bed > $intervene/${base}.bed
done

#Intersection of pooled peaks and replicated peaks was selected | and peaks in blacklisted area | 
# and peaks that intersect with large control peaks were removed to  yield final bed file
for i in $intervene/*.bed ; do
    echo $i;
    base=`basename -s .bed $i`
    bedtools intersect -a $MACS2Pooled/${base}_peaks.narrowPeak \
                       -b $intervene/${base}.bed -wa -u > $bedtools/${base}_rep.bed 
    bedtools intersect -a $bedtools/${base}_rep.bed  \
                       -b $Reference/hg38_blacklist.bed -v > $bedtools/${base}.bed 
done
wc -l $bedtools/*.bed 
# Because of rep1 ER_Veh too small
# cat $intervene/MCF7_Veh_ER/sets/10_MCF7_Veh_ER_rep1.bed $intervene/MCF7_Veh_ER/sets/11_MCF7_Veh_ER_rep1_MCF7_Veh_ER_rep2.bed > $intervene/MCF7_Veh_ER1.bed
# bedtools intersect -a $MACS2Pooled/MCF7_Veh_ER_peaks.narrowPeak -b $intervene/MCF7_Veh_ER1.bed -wa -u > $bedtools/MCF7_Veh_ER1_rep.bed 
# bedtools intersect -a $bedtools/MCF7_Veh_ER1_rep.bed -b $Reference/hg38_blacklist.bed -v > $bedtools/MCF7_Veh_ER1.bed 


# Visualization the overlap peaks
intervene venn -i $bedtools/MCF7_Veh_ER.bed $bedtools/MCF7_E2_ER.bed \
               --save-overlaps -o $intervene/ER_Veh_E2
intervene venn -i $bedtools/MCF7_Veh_FOXA1.bed $bedtools/MCF7_E2_FOXA1.bed \
               --save-overlaps -o $intervene/FOXA1_Veh_E2


intervene venn -i $bedtools/MCF7_Veh_ER.bed $bedtools/MCF7_Veh_FOXA1.bed \
               --save-overlaps -o $intervene/ER_FOXA1_Veh
intervene venn -i $bedtools/MCF7_E2_ER.bed $bedtools/MCF7_E2_FOXA1.bed \
               --save-overlaps -o $intervene/ER_FOXA1_E2







#Finding peak summits
# Many peaks have multiple summits.
for i in $MACS2Pooled/*_summits.bed ; do
    echo $i;
    base=`basename -s _summits.bed $i`
    bedtools intersect -a $bedtools/${base}.bed \
                       -b $MACS2Pooled/${base}_summits.bed \
                       -wa -wb > $bedtools/${base}_withsummits.bed;
done
wc -l $bedtools/*.bed

## Generate bed file for reference on deepTools analysis | Including ER_E2_Veh, ER_E2_only, ER_Veh_only

# Make bed files that have ER_only, ER_FOXA1, and FOXA1_only peaks.
bedtools intersect -a $bedtools/MCF7_Veh_ER1.bed -b $bedtools/MCF7_E2_ER.bed -u -wa > $bedtools/MCF7_Veh_E2_ER1.bed
bedtools intersect -a $bedtools/MCF7_Veh_ER1.bed -b $bedtools/MCF7_E2_ER.bed -v > $bedtools/MCF7_Veh_ER1_only.bed
bedtools intersect -a $bedtools/MCF7_E2_ER.bed -b $bedtools/MCF7_Veh_ER1.bed -v > $bedtools/MCF7_E2_ER1_only.bed

#And continues intersect with *_withsubmits.bed file
bedtools intersect -a $bedtools/MCF7_Veh_ER_withsummits.bed -b $bedtools/MCF7_Veh_ER_only.bed -wa > $bedtools/MCF7_Veh_ER_only_summits.bed
bedtools intersect -a $bedtools/MCF7_Veh_ER_withsummits.bed -b $bedtools/MCF7_Veh_E2_ER.bed -wa > $bedtools/MCF7_Veh_E2_ER_summits.bed
bedtools intersect -a $bedtools/MCF7_E2_ER_withsummits.bed -b $bedtools/MCF7_E2_ER_only.bed -wa > $bedtools/MCF7_E2_ER_only_summits.bed






wc -l $bedtools/MCF7_Veh_ER.bed
wc -l $bedtools/MCF7_E2_ER.bed


wc -l $bedtools/*_Veh_E2*.bed
wc -l $bedtools/*_only*.bed

mkdir -p $Results/deepTools/plotFingerprint \
         $Results/deepTools/multiBigwigSummary
plotFingerprint="$Results/deepTools/plotFingerprint"
multiBigwigSummary="$Results/deepTools/multiBigwigSummary"

mkdir -p $Results/deepTools/computeMatrix $Results/deepTools/plotHeatmap
matrixDir="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/deepTools/computeMatrix"
plotDir="/media/hkh/8TB/XUANTHANG/BreastCancer/ER_FOXA1_ChIPseq/Results/deepTools/plotHeatmap"

# WT and deltaC
computeMatrix reference-point --referencePoint center \
                              -R $bedtools/MCF7_Veh_ER_only_summits.bed \
                                 $bedtools/MCF7_Veh_E2_ER_summits.bed \
                                 $bedtools/MCF7_E2_ER_only_summits.bed \
                              -S $bedGraphToBigWig/ER_Veh.bw \
                                 $bedGraphToBigWig/ER_E2.bw \
                                 $bedGraphToBigWig/FOXA1_Veh.bw \
                                 $bedGraphToBigWig/FOXA1_E2.bw \
                              --missingDataAsZero -a 3000 -b 3000 --binSize 10 -p max \
                              -out $matrixDir/MCF7_E2_Veh.tab.gz

plotHeatmap -m $matrixDir/MCF7_E2_Veh.tab.gz -out $plotDir/MCF7_E2_Veh.pdf \
            --colorMap YlOrRd \
            --regionsLabel Loss, Common, Gain \
            --refPointLabel "Peak" 
# YlGnBu



#Make bed files with peaks containing location information of the summits using R
#file 




#Make combined peak file “.bed”.
cat $bedtools/MCF7_Veh_ER.bed $intervene/ER_Veh_E2/sets/01_MCF7_E2_ER.bed > $bedtools/MCF7_ER_Veh_E2.bed ## All peaks.

computeMatrix reference-point --referencePoint center -p max \
                              -R $bedtools/ER_Veh_E2_withsummits.bed \
                              -S $bedGraphToBigWig/ER_Veh.bw \
                                 $bedGraphToBigWig/ER_E2.bw \
                              --missingDataAsZero -a 1000 -b 1000 --binSize 10 \
                              -out $matrixDir/MCF7_E2_Veh_peaksummits.tab.gz

plotProfile -m $matrixDir/MCF7_E2_Veh_peaksummits.tab.gz \
            -out $plotDir/MCF7_E2_Veh_peaksummits.pdf --samplesLabel ER_Veh ER_E2 \
            --refPointLabel "ER Peak Summits" \
            --colors darkblue darkblue











