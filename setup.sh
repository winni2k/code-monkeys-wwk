#!/bin/bash

# download reference sequences
cd data
wget https://www.dropbox.com/s/swvm0zz0os3i1ji/hs37d5.chr20.fa.gz
wget https://www.dropbox.com/s/twed3rojdbhb1zg/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
cd ..

cd src
bash README.sh
cd ..
