

DATA_DIR=../../data
SAMP_DIR=$(DATA_DIR)/samples

REF_BASE=hs37d5.chr20
DICT=$(DATA_DIR)/$(REF_BASE).dict
REF=$(DATA_DIR)/$(REF_BASE).fa.gz
BIN=../../bin


SAMP=HG00096
REF_UZ=$(DATA_DIR)/$(REF_BASE).fa

all: varCalls/$(SAMP).raw.vcf

.PHONY: all

# create a bwa index of the human reference 
$(REF).bwt: $(REF)
	$(BIN)/bwa index $(REF)

# align HG00096 to GRCh37
bams/$(SAMP).bam: $(REF).bwt $(SAMP_DIR)/$(SAMP).fastq
	mkdir -p bams
	$(BIN)/bwa mem $(REF) $(SAMP_DIR)/$(SAMP).fastq -R '@RG\tID:foo\tSM:bar' | samtools view -b - > bams/$(SAMP).bam

#sort the bam
bams/$(SAMP).sorted.bam: bams/$(SAMP).bam
	$(BIN)/samtools sort bams/$(SAMP).bam -o bams/$(SAMP).sorted.bam -Obam -T tmp

#index the bam file
bams/$(SAMP).sorted.bam.bai: bams/$(SAMP).sorted.bam
	$(BIN)/samtools index bams/$(SAMP).sorted.bam

# uncompress the reference genome, as GATK only works on uncompressed fasta files
$(REF_UZ):$(REF)
	gzip -dc $(REF) > $(REF_UZ)

# index the uncompressed fasta file
$(REF_UZ).fai:$(REF_UZ)
	samtools faidx $(REF_UZ)

# create a dict of the fasta file
$(DICT):$(REF_UZ)
	java -jar $(BIN)/picard.jar CreateSequenceDictionary R=$(REF_UZ) O=$(DICT)

# call variants with the GATK
varCalls/$(SAMP).raw.vcf: bams/$(SAMP).sorted.bam $(REF_UZ) $(REF_UZ).fai $(DICT) bams/$(SAMP).sorted.bam.bai
	mkdir -p varCalls
	java -jar $(BIN)/GenomeAnalysisTK.jar \
	-R $(REF_UZ) \
	-T UnifiedGenotyper \
	-I bams/$(SAMP).sorted.bam \
	-o varCalls/$(SAMP).raw.vcf \
	-stand_call_conf 50.0 \
	-stand_emit_conf 10.0
