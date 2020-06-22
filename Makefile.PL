#Makefile.PL version 1.0

use 5.010001;
use ExtUtils::MakeMaker;
my @scripts = qw(htseq2multx.pl);

WriteMakefile
  (
   NAME            => 'htseq2multx',

   ($[ >= 5.005 ?  #Add these new keywords supported since 5.005

    (AUTHOR        => 'Robert W. Leach (rleach@princeton.edu)',
     ABSTRACT_FROM => "$scripts[0]") : ()),

   VERSION_FROM  => "$scripts[0]", # finds $VERSION
   PREREQ_PM     => {'IO::Pipe::Producer'     => 2.01,
                     'IO::Select'             => 1.21,
                     'IO::Uncompress::Gunzip' => 2.060,
                     'Getopt::Long'           => 2.38,
                     'File::Glob'             => 1.17,
                     'strict'                 => 1.07,
                     'warnings'               => 1.13,
                     'English'                => '1.06_01',
                     'Readonly'               => 2.05},
   EXE_FILES     => [@scripts]
  );
