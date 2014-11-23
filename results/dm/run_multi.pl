#!/bin/perl -w
use strict;

use DM;

use Getopt::Std;

my %args;
getopts( 'rj:', \%args );

my $dm = DM->new( dryRun => !defined $args{r}, numJobs => ( $args{j} || 1 ) );

my $DATA_DIR = '../../data';
my $SAMP_DIR = "$DATA_DIR/samples";

my $DICT = "$DATA_DIR/hs37d5.chr20.dict";
my $REF  = "$DATA_DIR/hs37d5.chr20.fa.gz";
my $BIN  = "../../bin";

# create a bwa index of the human reference
$dm->addRule( "$REF.bwt", $REF, "$BIN/bwa index $REF" );

# uncompress the reference genome, as GATK only works on uncompressed fasta files
my $REF_UZ = "$DATA_DIR/hs37d5.chr20.fa";
$dm->addRule( $REF_UZ, $REF, "gzip -dc $REF > $REF_UZ" );

# index the uncompressed fasta file
$dm->addRule( "$REF_UZ.fai", $REF_UZ, "samtools faidx $REF_UZ" );

# create a dict of the fasta file
$dm->addRule(
    $DICT,
    [ $REF_UZ, "$REF_UZ.fai" ],
    "rm -f $DICT"
      . " && java -jar $BIN/picard.jar CreateSequenceDictionary R=$REF_UZ O=$DICT"
);

for my $SAMP (qw/HG00096 HG00097 HG00099 HG00100 HG00101/) {

    # align sample to GRCh37
    $dm->addRule(
        "bams/$SAMP.bam",
        [ $REF, "$SAMP_DIR/$SAMP.fastq" ],
        "$BIN/bwa mem $REF $SAMP_DIR/$SAMP.fastq"
          . q/ -R '@RG\tID:foo\tSM:bar'/
          . " | samtools view -b - > bams/$SAMP.bam"
    );

    #sort the bam
    $dm->addRule( "bams/$SAMP.sorted.bam", "bams/$SAMP.bam",
"$BIN/samtools sort bams/$SAMP.bam -o bams/$SAMP.sorted.bam -Obam -T tmp"
    );

    #index the bam file
    $dm->addRule( "bams/$SAMP.sorted.bam.bai", "bams/$SAMP.sorted.bam",
        "$BIN/samtools index bams/$SAMP.sorted.bam" );

    # call variants with the GATK
    $dm->addRule(
        "varCalls/$SAMP.raw.vcf",
        [
            $REF_UZ, "$REF_UZ.fai",
            $DICT,   "bams/$SAMP.sorted.bam",
            "bams/$SAMP.sorted.bam.bai"
        ],
        "java -jar $BIN/GenomeAnalysisTK.jar"
          . " -R $REF_UZ"
          . " -T UnifiedGenotyper"
          . " -I bams/$SAMP.sorted.bam"
          . " -o varCalls/$SAMP.raw.vcf"
          . " -stand_call_conf 50.0"
          . " -stand_emit_conf 10.0"
    );

}

$dm->execute;
