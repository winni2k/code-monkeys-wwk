###########
# Sat Nov 22 19:46:27 GMT 2014
# downloaded bwa from http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download using browser
# unzip and compile bwa
tar xvjf bwa-0.7.10.tar.bz2 
cd bwa-0.7.10
make -j 2
install bwa ../../bin

# downloaded latest build of samtools from http://sourceforge.net/projects/samtools/
tar xvjf samtools-1.1.tar.bz2 
cd samtools-1.1/
make -j 2
install samtools ../../bin
