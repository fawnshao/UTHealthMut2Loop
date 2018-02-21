#!/bin/sh
# adapters=/home1/04935/shaojf/scratch/eCLIP_MCF-7/adapters.txt
# barcodes=/home1/04935/shaojf/scratch/eCLIP_MCF-7/primers.barcodes.txt
######## eCLIP adapters ########
# >RNA_A01
# AUUGCUUAGAUCGGAAGAGCGUCGUGUAG
# >RNA_B06
# ACAAGCCAGAUCGGAAGAGCGUCGUGUAG
# >RNA_C01
# AACUUGUAGAUCGGAAGAGCGUCGUGUAG
# >RNA_D08
# AGGACCAAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X1A
# AUAUAGGNNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X1B
# AAUAGCANNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X2A
# AAGUAUANNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RNA_X2B
# AGAAGAUNNNNNAGAUCGGAAGAGCGUCGUGUAG
# >RiL19
# AGAUCGGAAGAGCGUCGUG
# >AR17
# ACACGACGCTCTTCCGA
# >rand103T3
# NNNNNNNNNNAGATCGGAAGAGCACACGTCTG
# >PCR_F_D501
# AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
# >PCR_F_D502
# AATGATACGGCGACCACCGAGATCTACACATAGAGGCACACTCTTTCCCTACACGACGCTCTTCCGATCT
# >PCR_R_D701
# CAAGCAGAAGACGGCATACGAGATCGAGTAATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
# >PCR_R_D702
# CAAGCAGAAGACGGCATACGAGATTCTCCGGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC

# >Index_1_(i7)_Adapters
# CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
# >Index_2_(i5)_Adapters
# AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
# >Adapter_Trimming
# CTGTCTCTTATACACATCT

fastq_dir=/home1/04935/shaojf/scratch/eCLIP_MCF-7/fastqs

### for STAR
hg19file=/home1/04935/shaojf/scratch/star_index/hg19.fa
hg19indexdir=/home1/04935/shaojf/scratch/star_index/hg19.star
repfile=/home1/04935/shaojf/scratch/star_index/RepBase23.01.human.fa
repindexdir=/home1/04935/shaojf/scratch/star_index/RepBase23.01.human.star
### STAR index, run once, used in future
mkdir $hg19indexdir
mkdir $repindexdir
STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $hg19indexdir --genomeFastaFiles $hg19file
STAR --runMode genomeGenerate --runThreadN 68 --genomeDir $repindexdir --genomeFastaFiles $repfile

## step 1, FastQC
echo "step 1, FastQC ......"
mkdir fastqc_res
# fastqc --contaminants $barcodes --adapters $adapters --outdir fastqc_res $fastq_dir/*.fastq
fastqc --outdir fastqc_res $fastq_dir/*.fastq
### Demultiplexing
## don't need to do this
# for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
# do
# 	mkdir demux_paired_end_res
# 	demux_paired_end.py --fastq_1=$fastq_dir/${pre}_R1.fastq --fastq_2=$fastq_dir/${pre}_R1.fastq \
# 	--out_file_1=demux_paired_end_res/${pre}_R1.fastq --out_file_2=demux_paired_end_res/${pre}_R1.fastq
# done

## step 2, Cutadapt
echo "step 2, Cutadapt ......"
adapterTrimmed=fastqs_adapterTrim
mkdir $adapterTrimmed
# Cutadapt round 1: Takes output from demultiplexed files. Run to trim off both 5’ and 3’ adapters on both reads
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 --nextseq-trim=20 --cores=21 -A a0=NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g g1=CTTCCGATCTACAAGTT -g g2=CTTCCGATCTTGGTCCT -A a1=AACTTGTAGATCGGA -A a2=AGGACCAAGATCGGA -A a3=ACTTGTAGATCGGAA -A a4=GGACCAAGATCGGAA -A a5=CTTGTAGATCGGAAG -A a6=GACCAAGATCGGAAG -A a7=TTGTAGATCGGAAGA -A a8=ACCAAGATCGGAAGA -A a9=TGTAGATCGGAAGAG -A a10=CCAAGATCGGAAGAG -A a11=GTAGATCGGAAGAGC -A a12=CAAGATCGGAAGAGC -A a13=TAGATCGGAAGAGCG -A a14=AAGATCGGAAGAGCG -A a15=AGATCGGAAGAGCGT -A a16=GATCGGAAGAGCGTC -A a17=ATCGGAAGAGCGTCG -A a18=TCGGAAGAGCGTCGT -A a19=CGGAAGAGCGTCGTG -A a20=GGAAGAGCGTCGTGT -o $adapterTrimmed/${pre}_R1.adapterTrim.fastq -p $adapterTrimmed/${pre}_R2.adapterTrim.fastq $fastq_dir/${pre}_R1.fastq $fastq_dir/${pre}_R2.fastq > $adapterTrimmed/${pre}.fastq.adapterTrim.metrics &
done
wait
# Takes output from cutadapt round 1. Run to trim off the 3’ adapters on read 2, to control for double ligation events
# add "--nextseq-trim=20" for our experiment
# --cores=21 for fastq. 
# Make also sure that you have pigz (parallel gzip) installed if you use multiple cores and write to a .gz output file. Otherwise, compression of the output will be done in a single thread and therefore be the main bottleneck.
for pre in MCF7-CLIP-input MCF7-CLIP-RBM39-IPA MCF7-CLIP-RBM39-IPB
do
	cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 --cores=21 -A a1=AACTTGTAGATCGGA -A a2=AGGACCAAGATCGGA -A a3=ACTTGTAGATCGGAA -A a4=GGACCAAGATCGGAA -A a5=CTTGTAGATCGGAAG -A a6=GACCAAGATCGGAAG -A a7=TTGTAGATCGGAAGA -A a8=ACCAAGATCGGAAGA -A a9=TGTAGATCGGAAGAG -A a10=CCAAGATCGGAAGAG -A a11=GTAGATCGGAAGAGC -A a12=CAAGATCGGAAGAGC -A a13=TAGATCGGAAGAGCG -A a14=AAGATCGGAAGAGCG -A a15=AGATCGGAAGAGCGT -A a16=GATCGGAAGAGCGTC -A a17=ATCGGAAGAGCGTCG -A a18=TCGGAAGAGCGTCGT -A a19=CGGAAGAGCGTCGTG -A a20=GGAAGAGCGTCGTGT -o $adapterTrimmed/${pre}_R1.adapterTrim.round2.fastq -p $adapterTrimmed/${pre}_R2.adapterTrim.round2.fastq $adapterTrimmed/${pre}_R1.adapterTrim.fastq $adapterTrimmed/${pre}_R2.adapterTrim.fastq > $adapterTrimmed/${pre}.fastq.adapterTrim.round2.metrics &
done
wait

## step 3, STAR rmRep
echo "step 3, STAR rmRep ......"
adapterTrimmed=fastqs_adapterTrim
mkdir $adapterTrimmed
# STAR rmRep: Takes output from cutadapt round 2. Maps to human specific version of RepBase used to remove repetitive elements, helps control for spurious artifacts from rRNA (& other) repetitive reads.


STAR --runMode alignReads --runThreadN 16 --genomeDir
/path/to/RepBase_human_database_file --genomeLoad LoadAndRemove --
readFilesIn
/full/path/to/files/file_R1.C01.fastq.gz.adapterTrim.round2.fastq.gz
/full/path/to/files/file_R2.C01.fastq.gz.adapterTrim.round2.fastq.gz --
outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix
/full/path/to/files/file_R1.C01.fastq.gz.adapterTrim.round2.rep.bam --
outSAMattributes All --readFilesCommand zcat --outStd BAM_Unsorted --
outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped
Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType
EndToEnd >
/full/path/to/files/file_R1.C01.fastq.gz.adapterTrim.round2.rep.bam
