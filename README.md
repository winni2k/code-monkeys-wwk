Code Monkeys -- Using GNU make and DistributedMake (DM)
=======================================================

This is the git repository for the Code Monkeys session on GNU make
and the perl library DistributedMake (DM).  In this session we will
build a pipeline using a bash script, make file and a DM-based perl
script. Each of these methods have their advantages and
disadvantages. The goal of this session is to get a hands-on
understanding of what they are.

The pipeline will be borrowed from bioinformatics.  Specifically, we
will be aligning human Illumina sequencing reads to the human
reference genome, and processing the aligned reads with GATK.

Directory structure
--------

### data
Contains input data for the pipeline.  See `data/README.sh` for when
and where the files were downloaded.

### results
Contains the three working directories `results/bash`, `results/make`
and `results/dm`: one for each pipelining method.

### src
Contains all programs needed to perform the steps of the pipeline







