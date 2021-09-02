Fastq -> Trimming -> Alignment -> Remove duplicate -> MACS2 peakcalling (--keepdup all) -> Turn to R


======================================================================================================================
#1. Download rawdata, QC and Triming   
workDir="/media/hkh/8TB/XUANTHANG/BreastCancer/MCF7_E2_1h"
mkdir -p $workDir/RawData/sra $workDir/RawData/Fastq $workDir/RawData/FastQC $workDir/Results/Bowtie2 $workDir/Results/sambamba
export PATH=/media/hkh/8TB/XUANTHANG/TOOLS/sratoolkit.2.10.9-ubuntu64/bin:$PATH

RawData="$workDir/RawData"
Results="$workDir/Results"


# ChIPseq SINGLE-END
prefetch --option-file $RawData/SRP*.txt --output-directory $RawData/sra/
fastq-dump $RawData/sra/SRR*/SRR* --outdir $RawData/Fastq --gzip
fastqc -t 6 -o $RawData/FastQC --noextract -f fastq $RawData/Fastq/SRR*
multiqc -o $RawData -n MCF7_E2_1h_qc_report $RawData/FastQC/.


#Trimming by TrimGalore
mkdir $RawData/TrimGalore 
trim_galore -j 7 --illumina --fastqc --length 10 --gzip --output_dir $RawData/TrimGalore/ $RawData/Fastq/*.gz
multiqc $RawData/TrimGalore/. -o $RawData -n MCF7_E2_1h_Trim_report


=====================================================================================================================
#2 Alignment ChIPseq  (Single-end, Bowtie2 Using UCSC HUMAN hg38) 

# set up file names and grab base of filename for naming outputs
fq="$RawData/TrimGalore/*.fq.gz"
bowtie2Ref=/media/hkh/8TB/XUANTHANG/References/Refrence_Human_Genome/Homo_sapiens_UCSC_hg38/Sequence/Bowtie2Index/genome

#Bowtie2 alignment | samtools generated BAM file 
for i in $fq ; do
	echo $i;
	base=`basename -s _trimmed.fq.gz $i`
	bowtie2 -p 6 -x ${bowtie2Ref} -U $RawData/TrimGalore/${base}_trimmed.fq.gz 2> $Results/Bowtie2/${base}.log | \
	samtools view -@ 6 -h -S -b -o $Results/Bowtie2/${base}.bam
done


#Sort and index BAM file by genomic coordinates
for i in $Results/Bowtie2/*.bam ; do
	echo $i;
	base=`basename -s .bam $i`
	sambamba sort -t 6 -o $Results/Bowtie2/${base}_sorted.bam $Results/Bowtie2/${base}.bam 
done


#filtering Alignment file by sambamba or MarkDuplicate (Piscard - GATK tools) REMOVE_DUPLICATE=True
#The filter given are ‘not mapped’, ‘not duplicate’, and ‘[XS] ==null’, which are connected by ‘and’ operator.
# Filter out multi-mappers and duplicates

for i in $Results/Bowtie2/*_sorted.bam ; do
	echo $i;
	base=`basename -s _sorted.bam $i`
	sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" $Results/Bowtie2/${base}_sorted.bam > $Results/sambamba/${base}.bam
done

# Create indices for all the bam files for visualization and QC
for i in $Results/sambamba/*.bam ; do
	samtools index $i
done


#4 MACS2_Peakcalling 
mkdir -p $Results/MACS2Predict $Results/MACS2 $Results/MACS2log
Align="$Results/sambamba"
MACS2="$Results/MACS2"
log="$Results/MACS2log"

#4.1 Predict fragment length
for i in $Results/sambamba/*.bam ; do
	echo $i;
	base=`basename -s .bam $i`
	macs2 predictd -i $Results/sambamba/${base}.bam -g hs -m 5 20 2> $Results/MACS2Predict/${base}.txt
done




























#4.2 MACS2 callPeak
echo "****MACS2_Peakcalling in MCF7 celline****"
macs2 callpeak -t $Align/SRR2176971.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_Untreat_ER_rep1 -f BAM -B --nomodel --extsize 212 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_Untreat_ER_rep1.log
macs2 callpeak -t $Align/SRR2176972.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_Untreat_ER_rep2 -f BAM -B --nomodel --extsize 206 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_Untreat_ER_rep2.log
macs2 callpeak -t $Align/SRR2176973.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_E2_ER_rep1 -f BAM -B --nomodel --extsize 200 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_E2_ER_rep1.log
macs2 callpeak -t $Align/SRR2176974.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_E2_ER_rep2 -f BAM -B --nomodel --extsize 214 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_E2_ER_rep2.log
macs2 callpeak -t $Align/SRR2176975.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_Untreat_FOXA1_rep1 -f BAM -B --nomodel --extsize 181 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_Untreat_FOXA1_rep1.log
macs2 callpeak -t $Align/SRR2176976.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_Untreat_FOXA1_rep2 -f BAM -B --nomodel --extsize 159 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_Untreat_FOXA1_rep2.log
macs2 callpeak -t $Align/SRR2176979.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_E2_FOXA1_rep1 -f BAM -B --nomodel --extsize 174 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_E2_FOXA1_rep1.log
macs2 callpeak -t $Align/SRR2176980.bam -c $Align/SRR2176981.bam \
		-g hs -n MCF7_E2_FOXA1_rep2 -f BAM -B --nomodel --extsize 163 \
		--outdir $MACS2 --keep-dup all 2> $log/MCF7_E2_FOXA1_rep2.log

echo "****MACS2_Peakcalling in ZR751 celline****"
macs2 callpeak -t $Align/SRR2176986.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_Untreat_ER_rep1 -f BAM -B --nomodel --extsize 193 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_Untreat_ER_rep1.log
macs2 callpeak -t $Align/SRR2176987.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_Untreat_ER_rep2 -f BAM -B --nomodel --extsize 213 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_Untreat_ER_rep2.log
macs2 callpeak -t $Align/SRR2176988.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_E2_ER_rep1 -f BAM -B --nomodel --extsize 184 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_E2_ER_rep1.log
macs2 callpeak -t $Align/SRR2176989.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_E2_ER_rep2 -f BAM -B --nomodel --extsize 185 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_E2_ER_rep2.log
macs2 callpeak -t $Align/SRR2176990.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_Untreat_FOXA1_rep1 -f BAM -B --nomodel --extsize 211 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_Untreat_FOXA1_rep1.log
macs2 callpeak -t $Align/SRR2176991.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_Untreat_FOXA1_rep2 -f BAM -B --nomodel --extsize 218 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_Untreat_FOXA1_rep2.log
macs2 callpeak -t $Align/SRR2176994.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_E2_FOXA1_rep1 -f BAM -B --nomodel --extsize 239 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_E2_FOXA1_rep1.log
macs2 callpeak -t $Align/SRR2176995.bam -c $Align/SRR2176996.bam \
		-g hs -n ZR751_E2_FOXA1_rep2 -f BAM -B --nomodel --extsize 212 \
		--outdir $MACS2 --keep-dup all 2> $log/ZR751_E2_FOXA1_rep2.log

echo "****MACS2_Peakcalling in T47D celline****"
macs2 callpeak -t $Align/SRR2177001.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_Untreat_ER_rep1 -f BAM -B --nomodel --extsize 213 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_Untreat_ER_rep1.log
macs2 callpeak -t $Align/SRR2177002.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_Untreat_ER_rep2 -f BAM -B --nomodel --extsize 189 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_Untreat_ER_rep2.log
macs2 callpeak -t $Align/SRR2177003.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_E2_ER_rep1 -f BAM -B --nomodel --extsize 207 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_E2_ER_rep1.log
macs2 callpeak -t $Align/SRR2177004.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_E2_ER_rep2 -f BAM -B --nomodel --extsize 172 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_E2_ER_rep2.log
macs2 callpeak -t $Align/SRR2177005.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_Untreat_FOXA1_rep1 -f BAM -B --nomodel --extsize 213 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_Untreat_FOXA1_rep1.log
macs2 callpeak -t $Align/SRR2177006.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_Untreat_FOXA1_rep2 -f BAM -B --nomodel --extsize 162 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_Untreat_FOXA1_rep2.log
macs2 callpeak -t $Align/SRR2177009.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_E2_FOXA1_rep1 -f BAM -B --nomodel --extsize 209 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_E2_FOXA1_rep1.log
macs2 callpeak -t $Align/SRR2177010.bam -c $Align/SRR2177011.bam \
		-g hs -n T47D_E2_FOXA1_rep2 -f BAM -B --nomodel --extsize 150 \
		--outdir $MACS2 --keep-dup all 2> $log/T47D_E2_FOXA1_rep2.log


#5 Visualization
# Using bedGraphToBigWig tool (with treat.pileup.bdg) or bamCompare/bamCoverage from bedTools

mkdir -p $Results/deepTools $Results/bigWig/bigWiglogs

### Convert .bam files to .bigWig(.bw) files
for i in $Results/sambamba/*.bam ; do
	echo $i;
	base=`basename -s .bam $i`
	bamCoverage -p 6 --extendReads 150 --centerReads --binSize 10 --normalizeUsing CPM --smoothLength 60 \
	-b $Results/sambamba/${base}.bam \
	-o $Results/bigWig/${base}.bw 2> $Results/bigWig/bigWiglogs/${base}_bamCoverage.log
done

### 2.2. Run multiBigwigSummary to get average scores in genomic regions

# the average values for our test ENCODE ChIP-Seq datasets are computed for consecutive genome bins (default size: 10kb) by using the bins mode

multiBigwigSummary bins \
-b H3K4me1_chr22.bw H3K4me3_chr22.bw H3K27me3_chr22.bw H3K27ac_chr22.bw \
--labels H3K4me1 H3K4me3 H3K27me3 H3K27ac \
--region chr22 \
-out HMscores_per_bin.npz --outRawCounts HMscores_per_bin.tab
# Number of bins found: 5081








#6 countpeak and quality check ChIP signal
echo "****Count peaks called in each sample****"
wc –l $Results/MACS2/*.narrowPeak 

# sort peak by –log10(p-value)
echo "****Sort each of the narrowPeak files****" 
for i in $Results/MACS2/*.narrowPeak ; do
echo $i;
base=`basename -s .narrowPeak $i`
sort -k8,8nr $Results/MACS2/${base}.narrowPeak > $Results/MACS2/${base}_sorted.narrowPeak
done


#5 Handling replicates: Irreproducible Discovery Rate (IDR- for 2 replicates) or bedtools (more than 2 replicates)

mkdir $Results/idr

idr --samples $Results/MACS2/MCF7_Untreated_ER_rep*_sorted.narrowPeak \
--input-file-type narrowPeak --rank q.value -o $Results/idr/MCF7_Untreated_ER_idr -\
-plot -l $Results/idr/MCF7_Untreated_ER.idr.log












#Common peaks from each TF
wc -l Downstream_Analysis/idr/*_idr

# The total number of peaks that pass the threshold of IDR < 0.05
awk '{if($5 > 540) print $0}' Downstream_Analysis/idr/MCF7_Untreated_ER_idr | wc -l		
awk '{if($5 > 540) print $0}' Downstream_Analysis/idr/MCF7_E2_ER_idr | wc -l		


#7 R_code for Downstream analysis
cd /media/hkh/8TB/XUANTHANG/BreastCancer/FOXA1_ER_Cell2016
mkdir -p Downstream_Analysis/peaks Downstream_Analysis/reads
cp MACS2_Peakcalling/*.bed Downstream_Analysis/peaks
cp MACS2_Peakcalling/*sorted.narrowPeak Downstream_Analysis/peaks
mv Bowtie2Alignment/*.bam Downstream_Analysis/reads
mv Bowtie2Alignment/*.bai Downstream_Analysis/reads








