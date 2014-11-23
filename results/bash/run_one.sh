#!/bin/bash

DATA_DIR=../../data
SAMP_DIR=$DATA_DIR/samples

DICT=$DATA_DIR/hs37d5.chr20.dict
REF=$DATA_DIR/hs37d5.chr20.fa.gz
BIN=../../bin

# create a bwa index of the human reference 
$BIN/bwa index $REF

# align HG00096 to GRCh37
mkdir -p bams
SAMP=HG00096
$BIN/bwa mem $REF $SAMP_DIR/$SAMP.fastq -R '@RG\tID:foo\tSM:bar' | samtools view -b - > bams/$SAMP.bam

#sort the bam
$BIN/samtools sort bams/$SAMP.bam -o bams/$SAMP.sorted.bam -Obam -T tmp

#index the bam file
$BIN/samtools index bams/$SAMP.sorted.bam

# uncompress the reference genome, as GATK only works on uncompressed fasta files
REF_UZ=$DATA_DIR/hs37d5.chr20.fa
gzip -dc $REF > $REF_UZ

# index the uncompressed fasta file
samtools faidx $REF_UZ

# create a dict of the fasta file
java -jar $BIN/picard.jar CreateSequenceDictionary R=$REF_UZ O=$DICT

# call variants with the GATK
mkdir -p varCalls
java -jar $BIN/GenomeAnalysisTK.jar \
   -R $REF_UZ \
   -T UnifiedGenotyper \
   -I bams/$SAMP.sorted.bam \
   -o varCalls/$SAMP.raw.vcf \
   -stand_call_conf 50.0 \
   -stand_emit_conf 10.0

