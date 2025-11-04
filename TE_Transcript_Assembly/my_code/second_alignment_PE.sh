#! /bin/bash
#$1 input files name
#$2 STAR index position
#$3 threads

name=$1

STAR="/software/STAR-2.7.1a/bin/Linux_x86_64/STAR" # change this to STAR location

$STAR --runThreadN $3 \
  --genomeDir $2 \
  --readFilesIn  $name'_trimed_1.fastq' $name'_trimed_2.fastq'\
  --outFileNamePrefix  $name'_2nd_pass_'  \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 500 \
  --outSAMattributes NH HI NM MD XS AS \
  --genomeLoad LoadAndKeep \
  --alignIntronMax 1000000 \
  --outFilterType BySJout \
  --alignMatesGapMax 1000000 \
  --limitBAMsortRAM 30000000000
