#! /bin/bash
#$1 input files name
#$2 STAR index position
#$3 threads

name=$1
R1=$1'_R1.fastq.gz'
R2=$1'_R2.fastq.gz'

cutadapt="cutadapt" # change this to cutadapt location
STAR="STAR" # change this to STAR location

$cutadapt -j $3 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --quality-cutoff=15,10 --minimum-length=25  \
  -o $name'_trimed_1.fastq' -p  $name'_trimed_2.fastq'\
  $R1 $R2> $name'_cutadapt.trimlog'

$STAR --runThreadN $3 \
  --genomeDir $2 \
  --readFilesIn  $name'_trimed_1.fastq' $name'_trimed_2.fastq'\
  --outFileNamePrefix  $name'_1st_pass_'  \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 500 \
  --outSAMattributes NH HI NM MD XS AS \
  --genomeLoad LoadAndKeep \
  --alignIntronMax 1000000 \
  --outFilterType BySJout \
  --alignMatesGapMax 1000000 \
  --limitBAMsortRAM 30000000000
