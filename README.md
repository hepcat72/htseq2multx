# htseq2multx

Wrapper to [fastq-multx](https://github.com/brwnj/fastq-multx) to have it mimmic [barcode_splitter](https://bitbucket.org/princeton_genomics/barcode_splitter) and thus easily integrate with Princeton's [HTSEQ](http://htseq.princeton.edu) web interface.

## Description

The Princeton [HTSEQ](http://htseq.princeton.edu) software provides access to data generated by the [Genomics Core Facility](https://lsi.princeton.edu/facilities/sequencing-facility) in the [Lewis Sigler Institute for Integrative Genomics](https://lsi.princeton.edu).  It provides functionality for upstream analyses and data download/transfer.

Part of the supported functionality of [HTSEQ](http://htseq.princeton.edu) is the demultiplexing of pooled sequencing samples.  The initial implementation of this functionality utilized an in-house python tool called [barcode_splitter](https://bitbucket.org/princeton_genomics/barcode_splitter/) written in the interpreted scripting language Python.  The performance of barcode_splitter appears to be I/O bound, so as sample capacity grows and costs fall, it has proven to not scale well.

[fastq-multx](https://github.com/brwnj/fastq-multx) is a demultiplexing tool written in C and runs several orders of magnitude faster than [barcode_splitter](https://bitbucket.org/princeton_genomics/barcode_splitter).  Since [HTSEQ](http://htseq.princeton.edu) is highly integrated with [barcode_splitter](https://bitbucket.org/princeton_genomics/barcode_splitter), creating a wrapper to make [fastq-multx](https://github.com/brwnj/fastq-multx) look and act like [barcode_splitter](https://bitbucket.org/princeton_genomics/barcode_splitter) requires much less effort than altering [HTSEQ](http://htseq.princeton.edu) to use [fastq-multx](https://github.com/brwnj/fastq-multx).

## Install

This repo contains the details necessary to install the htseq2multx's module dependencies and script to both run the [fastq-multx](https://github.com/brwnj/fastq-multx) and convert its output.  It does not install [fastq-multx](https://github.com/brwnj/fastq-multx).  Refer to that [repo](https://github.com/brwnj/fastq-multx) for its installation instructions.

### System-wide install

    git clone https://github.com/hepcat72/htseq2multx.git
    cd htseq2multx
    perl Makefile.PL
    make
    sudo make install

### Local/user-account install

    git clone https://github.com/hepcat72/split-seq-processing.git
    cd split-seq-processing
    perl Makefile.PL INSTALL_BASE=~
    make
    make install

Where `~` in the perl command causes the module to go into ~/lib.  Replace `~` with wherever you wish to install.

### Caveats:

Module dependencies should install automatically, but if any do not and `perl Makefile.PL` issues a warning such as:

    Warning: prerequisite Readonly 2.05 not found.

you can try, for example:

    FTP_PASSIVE=1 PERL_MM_USE_DEFAULT=1 perl -MCPAN -e shell
    cpan> install Readonly

### Optional

    cd tests
    ./run_tests.sh
    make clean

You may need to either open a new terminal session or run the following before the executables will be recognized as being in your PATH:

- `hash -r` (bash)
- `rehash`  (tcsh)

## Usage

As an example, to split `forward_reads.fq` and `reverse_rease.fq` by sample, based on the sequences in `index.fq` and the barcodes defined in `barcodes.txt`, allowing 1 mismatch:

    htseq2multx.pl --mismatches 1 --bcfile barcodes.txt forward_reads.fq reverse_reads.fq index.fq --idxread 3

## LICENCE

See `LICENSE`