###########
# Sat Nov 22 18:48:22 GMT 2014
# downloaded the human reference genome chr20
wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz

# download 5 1000 genomes bam slices
mkdir -p samples
for i in HG00097 HG00099 HG00101; do
    ../bin/samtools view -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/$i/alignment/$i.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam 20:3000000-3100000 | ../bin/samtools bam2fq - > samples/$i.fastq &
done

for i in HG00096 HG00100; do
    ../bin/samtools view -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/$i/alignment/$i.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 20:3000000-3100000 | ../bin/samtools bam2fq - > samples/$i.fastq &
done



