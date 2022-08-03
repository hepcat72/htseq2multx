#!/usr/bin/env perl

use warnings;
use strict;

#This is to fool the executable version check that we are indeed fastq-multx
print("fastq-multx version 1.4\n");

if(scalar(@ARGV))
  {
    #Print an error that should cause a non-zero exit code
    print STDERR ("Returning error because of i/o error during file close\n");
  }

exit(0);
