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

Prerequisites
----
You will need a UNIX like environment. Furthermore, `Java JDK 7` or newer
and `GCC` are necessary to compile and run the programs that will be used in the
pipelines. For running the DistributedMake (DM) perl script located in
`results/dm` you will need `perl 5.8.8` or later as well as an
installed version of the perl library
[DM](https://github.com/wkretzsch/DM).  DM is still experimental and
has many perl library dependencies. As such, installing DM requires
advanced perl knowledge.  I highly recommend having a
working version of
[cpanm](http://search.cpan.org/~miyagawa/App-cpanminus-1.7016/lib/App/cpanminus.pm)
to install DM with. I believe DM version 0.13 is already installed on rescomp.

Getting started
----

Download this git repository

    git clone https://github.com/wkretzsch/code-monkeys-wwk.git
    cd code-monkeys-wwk

Run the setup script.  This script downloads a few large files I
deemed to large for git and compiles the programs in the `src`
directory

    ./setup.sh
    

Directory structure
--------

### data
Contains input data for the pipeline.  See `data/README.sh` for when
and where the files were downloaded.

### results
Contains the three working directories `results/bash`, `results/make`
and `results/dm`: one for each pipelining method.

### src
Contains all programs needed to perform the steps of the pipeline.

### bin
Compiled binaries are stored in `bin` for reference by the pipelines.





