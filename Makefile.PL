#Makefile.PL version 1.0

use strict;
use warnings;
use 5.010001;
use ExtUtils::MakeMaker;

my @scripts   = qw(htseq2multx.pl);
my $fqmx_link = "https://github.com/brwnj/fastq-multx";
my $fqmx_name = 'fastq-multx';

WriteMakefile
  (
   NAME            => 'htseq2multx',

   #Add these new keywords supported since 5.005
   ($] >= 5.005 ?

    (AUTHOR     => 'Robert W. Leach (rleach@princeton.edu)',
     ABSTRACT   => join('',
			('htseq2multx - Wrapper to fastq-multx to have it ',
			 'mimmic barcode_splitter and thus easily integrate ',
			 'with Princeton\'s HTSEQ web interface.'))) : ()),

   VERSION_FROM => "$scripts[0]",
   PREREQ_PM    => {'IO::Pipe::Producer'     => 2.01,
                    'IO::Select'             => 1.21,
                    'IO::Uncompress::Gunzip' => '2.060',
                    'Getopt::Long'           => 2.38,
                    'File::Glob'             => 1.17,
                    'strict'                 => 1.07,
                    'warnings'               => 1.13,
                    'English'                => 1.05,
                    'Readonly'               => 2.05},
   EXE_FILES    => [@scripts]
  );

my $fastq_multx = getExe($fqmx_name);
if(!defined($fastq_multx) || $fastq_multx eq '' ||
   !-e $fastq_multx || !-x $fastq_multx)
  {print STDERR ("WARNING: The $fqmx_name executable [$fastq_multx] appears ",
                 'to either not be in your path, not exist, not have execute ',
                 'permissions, or you have not created a symbolic link named ',
                 "'$fqmx_name' to the full name of the executable with ",
                 "version number.  If you have not installed $fqmx_name, ",
                 "you can find it here: $fqmx_link\n")}





sub getExe
  {
    my $command  = $_[0];
    my $sent_exe = $command;
    $sent_exe    =~ s/ .*//;
    my $exe      = '';

    if(eval("use File::Which;1;"))
      {
        $exe = which($sent_exe);
        #which from the module (at least in some versions) empties out valid
        #executables provided with absolute paths if that path is not in PATH,
        #in which case, this recovery is necessary
        if((!defined($exe) || $exe eq '') && -e $sent_exe && -x $sent_exe)
          {$exe = $sent_exe}
        elsif(!defined($exe))
          {$exe = ''}
      }
    else
      {
        print STDERR ("File::Which not found, switching to backup method.\n");
        $exe = `which $sent_exe`;
        chomp($exe);
        if($exe =~ /which: Command not found./)
          {
            print STDERR ("ERROR: System command 'which' does not appear to ",
                          'exist.  Please install the perl module File::Which.',
                          "\n");
            $exe = '';
          }
      }

    return($exe);
  }
