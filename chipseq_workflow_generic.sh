########################
# ChIP-seq data analysis scripts
# Author: Antonio Mora
# Created: Aug./2015
# Last modified: Mar./2016
########################
# NOTES:
# 1. This is a script written for a project comparing H2A.Z.1 and H2A.Z.2 binding profiles in mESC (paper to be submitted by Fosslie et al., 2016). It combines original code with useful lines of code adapted from multiple other sources.
# 2. The script has been written in order to be used in the Abel computer cluster of the University of Oslo (in case you don't have access to Abel, you can still adapt this script to your needs).
# 3. The script has been written for paired-end data, so it needs to be adapted for single-end data. It also assumes broad peaks, so it must be adapted for narrow peaks (macs2 statements). The reference genome used is mm9.
# 4. The script was applied to 32 samples relevant to the project (12 ours, 20 external). The following is a version of the script for two samples only.
########################

##################
# 1. Preparing directory and files (bash):
##################
WORKDIR=	# insert working directory here
cd $WORKDIR

# Copy "bedGraphToBigWig" to this folder:
wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz
gunzip v302.1.0.tar.gz
tar -xvf v302.1.0.tar
rm v302.1.0.tar
cp kentUtils-302.1.0/bin/linux.x86_64/bedGraphToBigWig .

# Download black-listed regions:
wget http://www.broadinstitute.org/~anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz
gunzip mm9-blacklist.bed.gz

# Unzip files:
gunzip *.fastq.gz

#################
# 2. QC analysis (bash):
#################
# Non-trimmed QC:
module load fastqc
fastqc *.fastq

# Trimming:
module load trim-galore

trim_galore --paired -o $WORKDIR INPUT_R1.fastq INPUT_R2.fastq
rm INPUT_R1.fastq INPUT_R2.fastq
trim_galore --paired -o $WORKDIR H2AZ_R1.fastq H2AZ_R2.fastq
rm H2AZ_R1.fastq H2AZ_R2.fastq

# Trimmed QC:
fastqc *.fq

##################
# 3. Alignment and Coverage (bash):
##################
module load bowtie2

# Indexing mm9:
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
gunzip chromFa.tar.gz
tar -xvf chromFa.tar
bowtie2-build chr1.fa,chr2.fa,chrX.fa,chr3.fa,chr4.fa,chr5.fa,chr7.fa,chr6.fa,chr8.fa,chr10.fa,chr14.fa,chr9.fa,chr11.fa,chr12.fa,chr13.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chrY.fa mm9_index

# Alignment Job (submit job to Slurm):
#==========================================
cat > antonio
#!/bin/bash
#SBATCH --job-name=antonio
#SBATCH --account= your-project-account
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1 --cpus-per-task=7

module load bowtie2
module load samtools
module load macs2

# INPUT:
bowtie2 -p 7 -x mm9_index -1 INPUT_R1_val_1.fq -2 INPUT_R2_val_2.fq -S INPUT.sam
rm INPUT_R1_val_1.fq INPUT_R2_val_2.fq
samtools view INPUT.sam -bS > INPUT.bam
samtools sort INPUT.bam INPUT_sorted
rm INPUT.bam

# H2AZ:
bowtie2 -p 7 -x mm9_index -1 H2AZ_R1_val_1.fq -2 H2AZ_R2_val_2.fq -S H2AZ.sam
rm H2AZ_R1_val_1.fq H2AZ_R2_val_2.fq
samtools view H2AZ.sam -bS > H2AZ.bam
samtools sort H2AZ.bam H2AZ_sorted
samtools mpileup H2AZ_sorted.bam  | perl -ne 'BEGIN{print "track type=wiggle_0 name=H2AZ description=H2AZ\n"};($c, $start, undef, $depth) = split; if ($c ne $lastC) { print "variableStep chrom=$c\n"; };$lastC=$c;next unless $. % 10 ==0;print "$start\t$depth\n" unless $depth<3;' > H2AZ.wig
macs2 callpeak -t H2AZ.sam -c INPUT.sam -n H2AZ --broad -g 1.87e9 -f SAM --broad-cutoff 0.01
rm H2AZ.bam
rm H2AZ.sam INPUT.sam
#==========================================
#(Ctrl+D)

sbatch antonio

# Alignment results per sample:
samtools flagstat INPUT_sorted.bam
samtools flagstat H2AZ_sorted.bam

# Number of peaks:
wc -l H2AZ_peaks.broadPeak

cat > cover
#!/bin/bash
#SBATCH --job-name=cover
#SBATCH --account= your-project-account
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1 --cpus-per-task=7

#Getting bedGraph files with genomeCoverageBed and sorting based on column 1, then column 2:
module load bedtools
genomeCoverageBed -ibam INPUT_sorted.bam -bg -g mm9.chrom.sizes | sort -k1,1 -k2,2n > INPUT.bedGraph
genomeCoverageBed -ibam H2AZ_sorted.bam -bg -g mm9.chrom.sizes | sort -k1,1 -k2,2n > H2AZ.bedGraph

sort -k1,1 -k2,2n INPUT.bedGraph > INPUT_sorted.bedGraph
sort -k1,1 -k2,2n H2AZ.bedGraph > H2AZ_sorted.bedGraph

#Making bigWig files from bedGraph:
./bedGraphToBigWig -unc INPUT_sorted.bedGraph mm9.chrom.sizes INPUT.bwig
./bedGraphToBigWig -unc H2AZ_sorted.bedGraph mm9.chrom.sizes H2AZ.bwig
#(Ctrl+D)

sbatch cover

########################
# 4. NGS-QC Quality control from bam files (galaxy):
########################
# Filezilla bam files sent to galaxy.ngs-qc.org

########################
# 5. Barplot of mapped vs unmapped reads (bash+R):
########################
# 5.1. Mapping results (bash):
samtools view -c -F 4 INPUT_sorted.bam	# s1_map
samtools view -c -f 4 INPUT_sorted.bam	# s1_unmap
samtools view -c -F 4 H2AZ_sorted.bam	# s2_map
samtools view -c -f 4 H2AZ_sorted.bam	# s2_unmap

# 5.2. Plots (R):
s1_map <-	# assign value of s1_map
s1_unmap <-	# assign value of s1_unmap
s2_map <-	# assign value of s2_map
s2_unmap <-	# assign value of s2_unmap

s1_tot <- s1_map + s1_unmap
s2_tot <- s2_map + s2_unmap

data <- matrix(c(s1_map/s1_tot, s1_unmap/s1_tot, s2_map/s2_tot, s2_unmap/s2_tot), nrow=2)
colnames(data) <- c("INPUT","H2AZ")

png("mapping_plot.png", height=8*300, width=8*300, res=300)
	barplot(data, main="Mapped to unmapped reads ratio per sample", xlab="H2A.Z sample", ylab="% reads", col=c("darkgreen","red"))
dev.off()

########################
# 6. Peak size analysis and Overlap analysis (R):
########################
# Get peak datasets:
setwd("your-working-directory")
H2AZ <- read.table("H2AZ_peaks.broadPeak", header=F, comment.char="", sep='\t', quote="")
H2AZ <- H2AZ[,1:4]
colnames(H2AZ) <- c("chr", "start", "end", "name")

# Get author's own comparison functions:
library(devtools)
install_github("antonio-mora/vennRanges")
library(vennRanges)

# Cleaning peak datasets:
mm9_black_list <- read.table("mm9-blacklist.bed", header=F, comment.char="", sep='\t', quote="")
colnames(mm9_black_list) <- c("chr", "start", "end")

H2AZ_clean <- range_comparison(H2AZ, mm9_black_list, "difference")

# Peak size analysis:
peak_stats <- function(bed) {
	num_peaks <- dim(bed)[1]
	size_peaks <- bed[,3] - bed[,2]
	res1 <- mean(size_peaks)
	res2 <- sd(size_peaks)
	res3 <- quantile(size_peaks, c(0.25,0.5,0.75,1))
	result <- list(num_peaks=num_peaks, mean_size=res1, sd_size=res2, quantiles=res3)
}

H2AZ_stats <- peak_stats(H2AZ_clean)

# Venn diagram analysis:
draw_2_venn(H2AZ_clean, other-file, "H2AZ_clean", "o-f", "red", "green")

########################
# 7. EXTRA. Additional analysis using tools that need to be installed (installation not explained here) (bash):
########################
# 7.1. Enrichment around Promoter and others (CEAS):
module load python2/2.7.10
module load R
export PYTHONPATH=$PYTHONPATH:$HOME/nobackup/software/ceas/lib/python2.7/site-packages/

cut -f1,2,3 H2AZ_peaks.broadPeak > H2AZ_peaks.bed
ceas -g mm9.refGene -b H2AZ_peaks.bed -w H2AZ.wig

# 7.2. RPKM plots (deepTools):
cat > deepTC
#!/bin/bash
#SBATCH --job-name=deepTC
#SBATCH --account= your-project-account
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1 --cpus-per-task=7

module load python2
module load samtools
export PATH=$PATH:/usit/abel/u1/antonimo/nobackup/software/deeptools/deepTools-2.2.3/bin

multiBamSummary bins --bamfiles S1.bam S2.bam S3.bam --binSize 1000 --labels S1 S2 S3 --outFileName inputs.npz

plotCorrelation --corData inputs.npz --plotFile correlation_plot_heat_inputs.pdf --corMethod spearman --whatToPlot heatmap --labels S1 S2 S3 --outFileCorMatrix correlation_matrix_heat_inputs.txt --colorMap YlOrRd --plotNumbers

plotCorrelation --corData inputs.npz --plotFile correlation_plot_scatter.pdf --corMethod spearman --whatToPlot scatterplot --labels S1 S2 S3 --outFileCorMatrix correlation_matrix_scatter.txt --colorMap Reds --plotNumbers --skipZeros --removeOutliers

plotPCA --corData inputs.npz --plotFile correlation_plot_PCA.pdf --labels S1 S2 S3 --outFileNameData PCA_data.txt


