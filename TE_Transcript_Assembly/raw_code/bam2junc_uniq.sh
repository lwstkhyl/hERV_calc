#!/bin/bash

# First parameter bam file
# Second parameter junction out

bamfile=$1
bedfile=$1.bed
juncfile=$2

samtools view -q 255 $bamfile | ~/transcript_assembly_pipeline_wanqing/filter_cs.py | ~/transcript_assembly_pipeline_wanqing/sam2bed.py > $bedfile
~/transcript_assembly_pipeline_wanqing/bed2junc.pl $bedfile $juncfile
