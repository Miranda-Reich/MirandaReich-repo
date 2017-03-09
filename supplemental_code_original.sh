#Supplemental code
#1. Fastqc quality assesment code

#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load gnu_parallel/201612222

###########  This is example code to move your files from our class shared to scratch, unzip them, concatenated the R1 and R2 files for each individual, and run fastqc the concatenated files 
## Recommended submission to ASC:
        # 10 cores
        # 2 hours
        # 20gb

######  remove any the targeted scratch directory and any files within
rm -r /scratch/aubcls06/DaphniaPunk
mkdir /scratch/aubcls06/DaphniaPunk

### Change directory to the scratch directory
cd /scratch/aubcls06/DaphniaPunk

#####   copy all .fastq.gz to  scratch directory in parallel using GNU parallel 
#  ls = make a (memory) list of the .fastq.gz files in that directory 
#  | in parallel using as many cores as possible (one job on each core) but no more jobs, cp the file in the list
# options: -- eta give the estimate time --dry-run to see if being parsed correctly
# You can do this seperately for the files you need to run using this code from before changing the specific file names

ls /home/aubcls06/class_shared/Exp1_DaphniaDataset/*.fastq.gz | time parallel -j+0 --eta 'cp {} .'

# unzip in parallel. List all the * .gz files and run on as many jobs as cores (-j) as and don't add any more files (0)
ls *fastq.gz |time parallel -j+0 'gunzip {}'

#### Create list of names:
# ls (list) contents of directory with fastq files, cut the names of the files at 
        #underscore characters and keep the first three chunks (i.e. fields; -f 1,2,3), 
        #sort names and keep only the unique ones (there will be duplicates of all 
        #file base names because of PE reads), then send the last 6 lines to a file 
        #called list with tail
                # 1 = C3
                # 2 = CCAGTT
                # 3 = L00
		# 4 = R1
                # 5 = 001

#### Make the list then use that list to Concatenate Forward Read (R1)files in parallel (Thanks Steven Sephic for this solution!)
ls | grep ".fastq" |cut -d "_" -f 1,2 | sort | uniq | parallel cat {1}_L00*_R1_*.fastq '>' {1}_All_R1.fastq ::: ${i}
##### Concatenate Reverse Reads (R2) files
ls | grep ".fastq" |cut -d "_" -f 1,2 | sort | uniq | parallel cat {1}_L00*_R2_*.fastq '>' {1}_All_R2.fastq ::: ${i}

##  Run fastqc on the All files in parallel
ls *_All_R1.fastq | time parallel -j+0 --eta 'fastqc {}'
ls *_All_R2.fastq | time parallel -j+0 --eta 'fastqc {}'

##Copy zipped files into my directory
cp *fastq.zip /home/aubcls06/class_shared/miranda_reich/Daphnia_data

#2. Trimming with Trimmomatic code

#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load trimmomatic/0.35
module load gnu_parallel/201612222

## Recommended submission:
        # 10 cores
        # 2 hours
        # 20gb


### Change directory to the scratch directory
cd /scratch/aubcls06/DaphniaPunk

##################  Now for the Cleaning ################################
# copy over the fasta file with the adapter file to use for screening
cp /home/aubcls06/class_shared/code/AdaptersToTrim_All.fa .

#### Create list of names:
# ls (list) contents of directory with fastq files, cut the names of the files at 
 #underscore characters and keep the first three chunks (i.e. fields; -f 1,2,3), 
        #sort names and keep only the unique ones (there will be duplicates of all 
        #file base names because of PE reads), then send the last 6 lines to a file 
        #called list with tail
                        # HS03_TTAGGC_L005_R1_001.fastq.gz 
                # 1 = C2
                # 2 = CCAGTT
                # 3 = All
                # 4 = R1
                # 5 = 001

ls | grep ".fastq" |cut -d "_" -f 1,2 | sort | uniq > list

### while loop to process through the names in the list
while read i
do

############ Trimmomatic #############
############  Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
#MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
#SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across
#requiredQuality: specifies the average quality required.
# -threads  is the option to define the number of threads (cores) to use. For this to be effective you need to request those cores at submission
#  ON HOPPER: trimmomatic-0.36

java -jar /opt/asn/apps/trimmomatic_0.35/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 10 -phred33 "$i"_All_R1.fastq "$i"_All_R2.fastq "$i"_All_R1_paired_threads.fastq "$i"_All_R1_unpaired_threads.fastq "$i"_All_R2_paired_threads.fastq "$i"_All_R2_unpaired_threads.fastq ILLUMINACLIP:AdaptersToTrim_All.fa:2:30:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

done<list

############### Now assess Quality again
#fastqc on the cleaned paired fastq files in parallel
ls *_R1_paired_threads.fastq | parallel -j+0  --eta 'fastqc {}'
ls *_R2_paired_threads.fastq | parallel -j+0  --eta 'fastqc {}'

##Copy zipped files into my directory
cp *_paired_threads_fastqc.zip /home/aubcls06/class_shared/miranda_reich/Daphnia_data

#3. Mapping and conversion from sam to bam using HiSat and Samtools

#!/bin/sh

#load modules
module load hisat2/2.0.5
module load samtools/1.3.1
module load gnu_parallel/201612222

### Change directory to the scratch directory
cd /scratch/aubcls06/Daphniapunk2

#Only way I could get this to run is when the index files are copied to DaphniaPunk scratch directory by hand.
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.1.ht2
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.2.ht2
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.3.ht2
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.4.ht2
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.5.ht2
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.6.ht2
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.7.ht2
# cp /home/aubcls06/class_shared/DaphniaPunk/Daphnia_pulex_INDEX3.8.ht2

#parallelize mapping
ls | grep "paired_threads.fastq" |cut -d "_" -f 1,2 | sort | uniq | time parallel -j+0 --eta hisat2 -p 10 --dta -x Daphnia_pulex_INDEX3 -1 {1}_All_R1_paired_threads.fastq {1}_All_R2_paired_threads.fastq -S {1}.sam ::: ${i}


while read i
do

samtools sort -@ 10 -o "$i".bam "$i".sam

done < list

#4. Assemble and merge transcripts together using StringTie

#!/bin/sh

#load modules
module load stringtie/1.3.3
module load gnu_parallel/201612222

#move into directory
cd /scratch/aubcls06/DaphniaPunk/

#copy over annotation file daphnia_genes2010_beta3.gtf before running script!!
#cp /home/aubcls06/class_shared/DaphniaPunk/daphnia_genes2010_beta3.gtf .

#Assemble transcripts for each sample
ls | grep ".bam" | cut -d "_" -f 1 | sort | uniq | time parallel -j+0 --eta stringtie -p 10 -G daphnia_genes2010_beta3.gtf -o {1}assembled.gtf -l {1} {1}*.bam ::: ${i}

#Get mergelist of assembled transcripts, letter 'a' added to transcript filename to allow grab of .gtf files ignoring index
ls *assembled.gtf | sort | uniq > mergelist.txt

#Merge transcripts from all samples
stringtie --merge -p 10 -G daphnia_genes2010_beta3.gtf -o stringtie_merged.gtf mergelist.txt

#Make appropriate directory heirarchy for ballgown
mkdir ballgown

#Estimate transcript abundances and create table counts for Ballgown
ls | grep ".bam" | cut -d "_" -f 1 | sort | uniq | time parallel -j+0 --eta stringtie -p 10 -G stringtie_merged.gtf -o ballgown/daph{1}/daph{1}.gtf -B -e {1}*.bam ::: ${i}

#Make tarball to bring over to R
tar -cvzf ballgown.tgz ballgown
cp *tgz /home/class_shared/aubcls06/miranda_reich/Daphnia_data

#5. Analysis differential gene expression using Ballgown in R

####Ballgown Differential Expression Analysis Basline####
#Code compiled from Pertea et. al 2016 Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown

#Load relevant R packages, downloaded previously
#Ballgown and genefilter DL from bioconductor, RSkittleBrewer from AlyssaFreeze github, and Devtools and dplyr installed using install.packages

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#Import the metadata summarizing the samples and treatments (must create this separately in a text editor;similar to colData file in garter snake practice)
pheno_data = read.csv("phenodata_Dp.csv")

#Import the countdata table made using the merge script.
# dataDir = the parent directory holding you StringTie output subdirectories/files
#samplePattern = pattern for ballgown to find in all your sample subdirectories
#daphnia_countdata = ballgown(dataDir = "CHANGEME", samplePattern = "CHANGEME", pData =pheno_data)
daphnia_countdata = ballgown(dataDir = "original_data", samplePattern = "daph", pData =pheno_data)

####Display data in tabular format####

#Filter low abundance genes out 
daphnia_countdata_filt = subset(daphnia_countdata, "rowVars(texpr(daphnia_countdata)) >1", genomesubset=TRUE)

#Identify transcripts that show sig difs between groups
results_transcripts = stattest(daphnia_countdata_filt, feature="transcript", covariate = "treatment", getFC = TRUE, meas="FPKM")

#Identify genes that show statistically sig diffs between groups
results_genes = stattest(daphnia_countdata_filt, feature="gene", covariate = "treatment", getFC = TRUE, meas="FPKM")

#Add gene names and gene IDs to the results_trascripts data frame
results_transcripts = data.frame(geneNames=ballgown::geneNames(daphnia_countdata_filt), geneIDs=ballgown::geneIDs(daphnia_countdata_filt), results_transcripts)

#Sort the results from smallest P value to the largest
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

#Write results to a csv file that can be shared and distributed
write.csv(results_transcripts, "daphnia_transcript_results.csv", row.names = FALSE)
write.csv(results_genes, "daphnia_genes_results.csv", row.names = FALSE)

#Identify transcripts and genes with a q value <0.05
sig.transcripts <- subset(results_transcripts,results_transcripts$qval<0.05)
sig.genes <- subset(results_genes,results_genes$qval<0.05)


####Display data in visual format####

#Make pretty plot colors
tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

#Show the distribution of gene abundances (measured as FPKM with ballgown) across samples, colored by treatment
fpkm = texpr(daphnia_countdata,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$treatment),las=2,ylab='log2(FPKM+1)')
#Plotting individual trascripts if you have one you want to target
ballgown::transcriptNames(daphnia_countdata)[12]
##  12
## "GTPBP6"
ballgown::geneNames(daphnia_countdata)[12]
##  12
## "GTPBP6
plot(fpkm[12,] ~ pheno_data$treatment, border=c(1,2),main=paste(ballgown::geneNames(daphnia_countdata)[12],' : ',ballgown::transcriptNames(daphnia_countdata)[12]),pch=19,xlab="Treatments",ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$treatment)),col=as.numeric(pheno_data$treatment))

#Plot structure and expression levels in a sample of all transcripts that share the same gene locus
plotTranscripts(ballgown::geneIDs(daphnia_countdata)[4.1], daphnia_countdata, main=c('Gene XIST in sample daphC3_CCGAAG'), sample=c('daphC3_CCGAAG'))

#Plot the average expression levels for all transcripts of a gene within different groups using the plotMeans function. Need to specify which gene to plot
plotMeans('MSTRG.7789', daphnia_countdata_filt, groupvar = "treatment",legend=FALSE)