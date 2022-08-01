#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long qw(:config no_auto_abbrev);
use IO::Pipe::Producer;
use IO::Select;
use English qw(-no_match_vars);
use Readonly;

our $VERSION     = '0.2';    #Version of this script
our $BCSVERSION  = '0.18.6'; #Emulated barcode_splitter version
our $FQMXVERSION = '1.4';    #Minimum required fastq-multx version

# CHANGE LOG
# 0.2 Added fatal error: Returning error because of i/o error during file close

#Exit codes (mimmicking barcode_splitter)
Readonly::Scalar my $SUCCESS         => 0;
Readonly::Scalar my $GENERIC_ERROR   => 1;
Readonly::Scalar my $OPEN_BCIN_ERROR => 5;
Readonly::Scalar my $OPEN_OUT_ERROR  => 6;
Readonly::Scalar my $OUTFILES_EXIST  => 9;
Readonly::Scalar my $FQIN_OPEN_ERROR => 10;
Readonly::Scalar my $DUPE_BC_ROWS    => 12;

#Number of fastq lines at which to give up trying to set the defline separator
Readonly::Scalar my $FORMAT_DETECT_LINE_MAX => 100;

my $help            = 0;
my $version         = 0;
my $idxread         = [];
my $mismatches      = 1;
my $barcodes_at_end = 0;
my $galaxy          = 0;
my $sanitize        = 1;
my $verbose         = 0;
my $gzipout         = 0;
my $split_all       = 0;
my $fastq_files     = [];
my $index_hash      = {};
my $bchash          = {};
my $samplehash      = {};
my $fastq_multx     = 'fastq-multx';
my $debug           = 0;
my $check_deflines  = 0;
my $tmp_files       = [];
my($bcfile,$bcfilefqmx,$prefix,$suffix,$format,$gzipin,$deflinesep);

main();

sub main
  {
    processOptions();

    my $command = getFastqMultxCommand();

    if($debug)
      {print STDERR ("#$command\n")}

    runFastqMultx($command);
  }






##
## Methods
##

sub processOptions
  {
    ##
    ## Parse the command line and handle help/version options
    ##

    my $OptionHash = {
                      'version:+'       => \$version,
                      'h|help:+'        => \$help,
                      'bcfile=s'        => \$bcfile,
                      'idxread=i{1,2}'  => $idxread,
                      'mismatches=i'    => \$mismatches,
                      'barcodes_at_end' => \$barcodes_at_end,
                      'prefix=s'        => \$prefix,
                      'suffix=s'        => \$suffix,
                      'galaxy'          => \$galaxy,
                      'sanitize!'       => \$sanitize,
                      'v|verbose:+'     => \$verbose,
                      'gzipout'         => \$gzipout,
                      'split_all'       => \$split_all,
                      'format=s'        => \$format,
                      'gzipin'          => \$gzipin,
                      '<>'              => sub {push(@$fastq_files,$ARG[0])},
                      'fastq-multx=s'   => \$fastq_multx,
                      'debug:+'         => \$debug,
                      'check-deflines'  => \$check_deflines
                     };

    if(!GetOptions(%$OptionHash))
      {
        print STDERR ("Unable to parse command line.\n");
        exit($GENERIC_ERROR);
      }

    ##
    ## Handle the version option
    ##

    #Emulate barcode_splitter
    if($version == 1)
      {print("$BCSVERSION\n")}
    #Provide all version information
    elsif($version)
      {
        my($multx_version,$is_fqmx,$problem) = getMultxVersion($fastq_multx);
        my $multx_used = ($is_fqmx ? join('.',@$multx_version) : 'unknown');
        print("htseq2multx, version $VERSION\n",
              "barcode_splitter, emulated version $BCSVERSION\n",
              "fastq-multx, minimum required version $FQMXVERSION\n",
              "$fastq_multx, version used $multx_used\n");
      }
    if($version)
      {exit($SUCCESS)}

    #Handle the help option
    if($help == 1)
      {help() && exit($SUCCESS)}
    elsif($help)
      {help2() && exit($SUCCESS)}

    ##
    ## Ensure basically valid values for required options are supplied
    ##

    if(scalar(@$fastq_files) == 0)
      {
        help(1);
        print STDERR ("barcode_splitter: error: the following arguments are ",
                      "required: FILE\n");
        exit($GENERIC_ERROR);
      }

    if(scalar(@$idxread) == 0)
      {
        help(1);
        print STDERR ("barcode_splitter: error: Sequence files and at least ",
                      "one number indicating the indexed file(s) (--idxread) ",
                      "is required\n");
        exit($GENERIC_ERROR);
      }

    if(!defined($bcfile))
      {
        help(1);
        print STDERR ('barcode_splitter: error: Must specify a barcodes file ',
                      'with "--bcfile" option',"\n");
        exit($GENERIC_ERROR);
      }
    elsif(!-e $bcfile)
      {
        #barcode_splitter does not print: help(1);
        print STDERR ("ERROR: Unable to open barcode file: [Errno 2] No such ",
                      "file or directory: '$bcfile'\n");
        exit($OPEN_BCIN_ERROR);
      }

    if(scalar(grep {$ARG < 1 || $ARG > scalar(@$fastq_files)} @$idxread))
      {
        help(1);
        print STDERR ('barcode_splitter: error: Invalid index read number ',
                      '("--idxread"), must be between 1 and 1 (the number of ',
                      "supplied sequence files)\n");
        exit($GENERIC_ERROR);
      }

    ##
    ## Validate the fastq-multx executable
    ##

    #Validate the executable (partially for security purposes, since this
    #wrapper may be initiated via a web tool)
    my ($multx_version,$is_fqmx,$problem) = getMultxVersion($fastq_multx);
    my $supported_version                 = [split(/\./,$FQMXVERSION)];
    if($problem =~ /./)
      {
        help2();
        print('ERROR: There is a problem with the fastq-multx executable: ',
              "[$problem].  Use the --fastq-multx option described above to ",
              "supply fastq-multx if it is missing from your PATH.\n");
        exit($GENERIC_ERROR);
      }
    elsif(!$is_fqmx)
      {
        help2();
        print('ERROR: The fastq-multx executable appears to not be fastq-',
              'multx according to the first 2 lines of its usage.  Use the ',
              "--fastq-multx option described above to supply fastq-multx.\n");
        exit($GENERIC_ERROR);
      }
    elsif(fqmxVersionInsufficient($multx_version,$supported_version))
      {
        help2();
        if(scalar(@$multx_version) < 1)
          {print STDERR ('Unable to determine fastq-multx version of ',
                         "executable [$fastq_multx].  Unable to proceed.\n")}
        else
          {print STDERR ('Version [',join('.',@$multx_version),
                         "] of the fastq-multx executable [$fastq_multx] is ",
                         "unsupported.  Only version [$FQMXVERSION] or higher ",
                         'is supported.  Use the --fastq-multx option ',
                         "described above to supply fastq-multx.\n")}
        exit($GENERIC_ERROR);
      }

    #Initialize the data structures/hashes used, e.g. index_hash, bchash,
    #samplehash, and read in & convert the barcode file
    my $order = 1;
    foreach my $fqf_num (@$idxread)
      {$index_hash->{$fqf_num - 1} = $order++}
    ($bcfilefqmx,$bchash,$samplehash) =
      bcs2fqmxBCFile($bcfile,scalar(@$idxread));

    ##
    ## Process the prefix & suffix values
    ##

    if(!defined($prefix))
      {$prefix = ''}

    if(defined($prefix) && $galaxy)
      {$prefix =~ s/_/-/g}

    if(defined($suffix) && $galaxy)
      {$suffix =~ s/_/-/g}

    if(!defined($gzipin))
      {$gzipin = ($fastq_files->[0] =~ /\.gz/)}

    if(!defined($gzipout))
      {
        if(defined($suffix))
          {$gzipout = ($suffix =~ /\.gz$/)}
        else
          {$gzipout = $gzipin}
      }
    elsif($gzipout)
      {
        if(defined($suffix) && $suffix !~ /\.gz$/)
          {$suffix .= '.gz'}
        elsif(!defined($suffix))
          {$suffix = '.fastq.gz'}
      }

    if(!defined($suffix))
      {$suffix = ($gzipout ? '.gz' : '.fastq')}

    ##
    ## Process the fastq file inputs
    ##

    #fastq-multx doesn't deal with gzipped input that doesn't have .gz in the
    #name, so those files must be written to an actual file
    if($gzipin)
      {$tmp_files = processGZInfiles($fastq_files,$prefix)}

    if($check_deflines)
      {
        $deflinesep = getDeflineSep($fastq_files);
        if($deflinesep eq '')
          {
            if(scalar(@$fastq_files) > 1)
              {print STDERR ('WARNING: Cannot identify fastq record type.  ',
                             "Unable to check that deflines match.\n")}
            $check_deflines = 0;
          }
      }

    my $missing_infiles = [grep {!-e $ARG} @$fastq_files];
    if(scalar(@$missing_infiles))
      {
        print STDERR ('ERROR: ',scalar(@$missing_infiles),' input files do ',
                      'not exist: ',join(',',@$missing_infiles),"\n");
        exit($FQIN_OPEN_ERROR);
      }

    ##
    ## Check the output files
    ##

    if(!$split_all && scalar(@$fastq_files) == scalar(@$idxread))
      {$split_all = 1}

    if(outfilesExist($bchash,$fastq_files,$split_all))
      {exit($OUTFILES_EXIST)}

    if(isOutdirMissing($prefix))
      {
        print STDERR ("ERROR: Directory in prefix [$prefix] does not exist.\n");
        exit($GENERIC_ERROR);
      }
  }

sub help
  {
    my $short_only = $ARG[0] || 0;

    my $short_msg  = << 'END_SHORT';
usage: htseq2multx.pl   [-h] [--version] [--bcfile FILE]
                        [--idxread READNUM [READNUM ...]]
                        [--mismatches MISMATCHES] [--barcodes_at_end]
                        [--prefix PREFIX] [--suffix SUFFIX] [--galaxy]
                        [--nosanitize] [-v] [--gzipout] [--split_all]
                        [--format FORMAT] [--gzipin]
                        FILE [FILE ...]
END_SHORT

    my $long_msg = << 'END_LONG';

Split one or more fastq files based on barcode sequence.

optional arguments:
  -h, --help            show this help message and exit*
  --version             show program's version number and exit*

Barcodes:
  --bcfile FILE         REQUIRED: Tab delimited file: "Sample_ID <tab>
                        Barcode_Sequence" Multiple barcode columns with
                        different barcode lengths allowed, but all barcodes in
                        each inidividual column must be the same length.*
  --idxread READNUM [READNUM ...]
                        REQUIRED: Indicate in which read file(s) to search for
                        the corresponding column of barcode sequences, e.g. if
                        the first column of barcodes is in the second sequence
                        read file and the second column's barcodes are in the
                        third sequence read file, you'd supply `--idxread 2 3`
  --mismatches MISMATCHES
                        Number of mismatches allowed in barcode matching
  --barcodes_at_end     Barcodes are at the end of the index read (default is
                        at the beginning)

Output Options:
  --prefix PREFIX       Prefix for output files
  --suffix SUFFIX       Suffix for output files (default based on --format)
  --galaxy              Produce "Galaxy safe" filenames by removing
                        underscores (default: False)
  --nosanitize          Do not produce "safe" filenames by replacing unusual
                        characters in the supplied prefix and sample IDs with
                        underscores. (default: False)
  -v, --verbose         verbose output
  --gzipout             Output files in compressed gzip format (default is
                        uncompressed)
  --split_all           Split all input files, including index read files (by
                        default, index read files are not split unless all
                        read files are index files)

Input format:
  --format FORMAT       Specify format for sequence files (fasta or fastq)
  --gzipin              Assume input files are in gzip format, despite file
                        extension (default is auto based on input file
                        extension)

Sequence Inputs:
  FILE                  A series of 1 or more [optionally zipped] fastq files.

* Supply --help twice to see differences between barcode_splitter and this tool.
END_LONG

    print($short_msg);
    unless($short_only)
      {print($long_msg)}
  }

sub help2
  {
    my $msg = << 'END_HELP2';

The following extra and modified options apply specifically to the htseq2multx.pl wrapper (as opposed to barcode_splitter).

arguments with behavioral differences:
  --bcfile FILE         barcode_splitter requires all barcodes in a single
                        column to be the same length, however fastq-multx
                        supports variable length barcodes, thus the length
                        consistency check is not performed and no errors about
                        different barcode lengths are issued.

optional arguments:
  -h, --help            Supply twice to show this message and exit.
  --version             Supply once to get the barcode_splitter version that
                        this wrapper emulates.  Supply twice to get a more
                        verbose version of the barcode_splitter version that
                        this wrapper emulates and the version of this wrapper.

configuration arguments:
  --fastq-multx PATH    Path to the fastq-multx installation.  Default is to use
                        the PATH in the user's environment.
  --debug               Supply this flag to output the fastq-multx command on
                        standard error and not delete the temporary barcode
                        file.
END_HELP2

    print($msg);
  }

sub getFastqMultxCommand
  {
    # -x          Don't trim barcodes off before writing out destination
    # -d N        Require a minimum distance of N between the best and next best (2)
    # -B BCFIL    Use barcodes from BCFIL, no determination step, codes in <read1.fq>
    # -m N        Allow N mismatches in union of all indexes, unless -M is supplied. (1)
    # -M M        Allow N,M mismatches in indexes 1,2 respectively (see -m N) (not used)
    my $command =
      "$fastq_multx -x -d 1 -B '$bcfilefqmx' -m $mismatches -M $mismatches";

    $command .= ($barcodes_at_end ? ' -e' : ' -b');
    $command .= ($check_deflines  ? " -v '$deflinesep'" : '');

    #Add the input files and output file templates
    my $templates_str = '';
    foreach my $fqf_index (sort {indexesFirst($a,$b)}
                           0..$#{$fastq_files})
      {
        $command       .= " '$fastq_files->[$fqf_index]'";
        $templates_str .= " -o '" . getOutfileTemplate($fqf_index) . "'";
      }
    $command .= $templates_str;

    return($command);
  }

#Runs the command and prints stdout/stderr.  Exits as barcode_splitter would.
sub runFastqMultx
  {
    my $command = $ARG[0];

    my $producer = IO::Pipe::Producer->new();
    my($stdout_handle,$stderr_handle) = $producer->getSystemProducer($command);

    my $sel = IO::Select->new();
    $sel->add($stdout_handle,$stderr_handle);
    my $stdout = '';
    my $make_error_fatal = 0;
    while(my @fhs = $sel->can_read())
      {
        foreach my $fh (@fhs)
          {
            my $line = <$fh>;
            unless(defined($line))
              {
                $sel->remove($fh);
                next;
              }
            if($fh == $stdout_handle)
              {$stdout .= $line}
            elsif($fh == $stderr_handle)
              {
                #Echo STDERR and exit non-0 if error is fatal
                if(processSTDERR($line))
                  {$make_error_fatal = 1}
              }
          }
      }

    if($CHILD_ERROR || $make_error_fatal)
      {
        print STDERR ("ERROR: fastq-multx command failed",
                      ($CHILD_ERROR ? ' with a non-zero exit code ' .
                       "[$CHILD_ERROR]" .
                       ($OS_ERROR =~ /./ ? " and error: $OS_ERROR" : '.') :
                       ':'),"\n\t$command\n");
        unless($stdout =~ /\t/)
          {exit($GENERIC_ERROR)}
      }

    print(fqmx2bcsStdout($stdout,scalar(@$idxread)));
  }

#Since fastq-multx does not handle gzipped input from process substitution, such
#input is detected here and written to temporary files.  Assumes $gzip is true.
sub processGZInfiles
  {
    my $fastq_files = $ARG[0];
    my $prefix      = $ARG[1];
    if(scalar(grep {$ARG !~ /\.gz$/i} @$fastq_files) == 0)
      {return()}

    my $tmp_files = [];

    #Create a temporary set of files
    my $time = time();
    my $dir = '.';
    if(defined($prefix) && $prefix =~ m%(.*/)%)
      {$dir = $1}
    my $template = "$dir/${time}-\%.fq.gz";

    for(my $i = 0;$i <= $#{$fastq_files};$i++)
      {
        if($fastq_files->[$i] =~ /\.gz$/i)
          {next}
        my $file = $template;
        $file =~ s/%/$i/;
        push(@$tmp_files,$file);
        writeGZProcSubFile($fastq_files->[$i],$file);
        $fastq_files->[$i] = $file;
      }

    return(wantarray ? @$tmp_files : $tmp_files);
  }

sub writeGZProcSubFile
  {
    my $inprocsub = $ARG[0];
    my $outfile   = $ARG[1];

    if(!open(IN,'<',$inprocsub))
      {
        print STDERR ('ERROR: Unable to open input stream from process ',
                      "substitution [$inprocsub].\n$OS_ERROR\n");
        exit($FQIN_OPEN_ERROR);
      }
    if(!defined(binmode(IN)))
      {
        print STDERR ('ERROR: Unable to put input stream from process ',
                      "substitution [$inprocsub] in binary mode.\n$OS_ERROR\n");
        exit($FQIN_OPEN_ERROR);
      }
    if(!open(OUT,'>',$outfile))
      {
        print STDERR ("ERROR: Unable to open temporary output file [$outfile] ",
                      "for process substitution [$inprocsub].\n$OS_ERROR\n");
        exit($OPEN_OUT_ERROR);
      }
    if(!defined(binmode(OUT)))
      {
        print STDERR ("ERROR: Unable to put temporary output file [$outfile] ",
                      "for process substitution [$inprocsub] in binary mode.\n",
                      "$OS_ERROR\n");
        exit($OPEN_OUT_ERROR);
      }

    print OUT (<IN>);

    if(!close(OUT))
      {
        print STDERR ("ERROR: Unable to close temporary output file [$outfile]",
                      " for process substitution [$inprocsub].\n$OS_ERROR\n");
        exit($GENERIC_ERROR);
      }
    if(!close(IN))
      {
        print STDERR ('ERROR: Unable to close input stream from process ',
                      "substitution [$inprocsub].\n$OS_ERROR\n");
        exit($FQIN_OPEN_ERROR);
      }
  }

sub isOutdirMissing
  {
    my $prefix  = $ARG[0];
    my $missing = 0;
    if($prefix =~ m%(.*/)%)
      {
        my $dir = $1;
        if(!-e $dir)
          {$missing = 1}
      }
    return($missing);
  }

sub outfilesExist
  {
    my $hash     = $ARG[0];
    my $files    = $ARG[1];
    my $existing = 0;

    foreach my $fqf_index (0..$#{$files})
      {
        my $template = getOutfileTemplate($fqf_index);
        if($template =~ m%n/a%i)
          {next}
        foreach my $sample (keys(%$hash),'unmatched')
          {
            my $outfile = $template;
            $outfile =~ s/\%/$sample/;
            if(-e $outfile)
              {
                print STDERR "ERROR: File exists: $outfile.\n";
                $existing = 1;
              }
          }
      }

    return($existing);
  }

sub getMultxVersion
  {
    my $version = [];
    my $error   = '';
    my $is_fqmx = 1;

    my $output  = `($fastq_multx 2>&1 ) | head -n 2`;
    if($CHILD_ERROR)
      {
        $error = "fastq-multx failed with a non-zero exit code [$CHILD_ERROR]" .
          ($OS_ERROR =~ /./ ? " and error: $OS_ERROR" : '.') . "\n";
        return($version,$is_fqmx,$error);
      }
    elsif($output !~ /fastq-multx/)
      {
        $is_fqmx = 0;
        return($version,$is_fqmx,$error);
      }
    elsif($output =~ /version\D*([\d\.]+)/i)
      {
        my $version_string = $1;
        push(@$version,split(/\./,$version_string));
      }
    return($version,$is_fqmx,$error);
  }

sub fqmxVersionInsufficient
  {
    my $used = $ARG[0];
    my $reqd = $ARG[1];

    my $lim =
      (scalar(@$used) < scalar(@$reqd) ? scalar(@$used) : scalar(@$reqd)) - 1;

    if($lim == 0)
      {return(1)}

    my $insufficient = 0;
    my $num_equal    = 0;
    foreach my $i (0..$lim)
      {
        if($used->[$i] < $reqd->[$i])
          {
            $insufficient = 1;
            last;
          }
        elsif($used->[$i] > $reqd->[$i])
          {last}
        else
          {$num_equal++}
      }

    #Add 1 because num_equal is 1-based and lim is 0-based
    $lim++;

    #If the required version matches, but the version number of the requirement
    #is longer than the numbers that were checked
    if(!$insufficient && $lim == $num_equal && scalar(@$reqd) > $num_equal)
      {$insufficient = 1}

    return($insufficient);
  }

#Convert a barcode file from barcode_splitter's accepted format to fastq-multx.
#  It additionally checks the sample names & barcodes and sanitizes the sample
#  names and it creates 2 hashes: $bchash and $samplehash.  It returns the
#  converted barcode file and the 2 hashes.  Example:
#
#($converted_file,$bchash,$samplehash) =
#  bcs2fqmxBCFile($barcode_splitter_barcode_file,$number_of_bc_columns)
#
#bchash is used to add barcodes to the output stats on STDOUT & looks like this:
#  bchash->{$sanitized_sample_name} = tab-delimited string of barcodes
#samplehash maps the sanitized sample name to their unsanitized version, e.g.:
#  samplehash->{$sanitized_sample_name} = $original_sample_name
sub bcs2fqmxBCFile
  {
    my $bcs_file    = $ARG[0];
    my $num_indexes = $ARG[1];
    my $fqmx_file   = $bcs_file . '.fqmx';
    my $bchash      = {};
    my $samplehash  = {unmatched => 'unmatched'};
    my $uniq_check  = {};
    my $dupes       = 0;

    unless(open(BCSBCF,$bcs_file))
      {print STDERR ("ERROR: Unable to open barcode file: $OS_ERROR") &&
         exit($OPEN_BCIN_ERROR)}

    unless(open(FQMXBCF,'>',$fqmx_file))
      {print STDERR ("ERROR: Unable to open fqmx barcode file: $OS_ERROR") &&
         exit($GENERIC_ERROR)}

    my $ln     = 0;
    while(<BCSBCF>)
      {
        $ln++;
        chomp;
        if(/^\s*#/ || /^\s*$/)
          {next}
        my @cols = split(/\t/,$ARG);
        if(scalar(@cols) >= ($num_indexes + 1))
          {
            #Keep track of sample names so they can be changed back
            #Barcode splitter only sanitizes/makes-galaxy-friendly the file
            #names, not sample names
            my $origsn = $cols[0];

            if($sanitize)
              {$cols[0] =~ s/[^A-Za-z0-9\-\+\.<>\@#\%^=_~]+/_/g}
            #These minimal changes are required for fastq-multx to not error out
            else
              {$cols[0] =~ s/[\/\s]/_/g}
            if($galaxy)
              {$cols[0] =~ s/_/-/g}

            $samplehash->{$cols[0]} = $origsn;

            my $bc_str = join('-',@cols[1..$num_indexes]);
            if(exists($uniq_check->{$bc_str}))
              {$dupes = 1}
            push(@{$uniq_check->{$bc_str}},$cols[0]);

            print FQMXBCF ("$cols[0]\t$bc_str\n");

            if(exists($bchash->{$cols[0]}))
              {
                print STDERR ("Sample IDs are not unique.  Note, this could ",
                              "be due to either being in --galaxy mode or ",
                              "from sample ID sanitization.  If the sample ",
                              "IDs in your barcodes file appear unique, ",
                              "either replace UTF8 characters with ASCII ",
                              "characters or try supplying --nosantize.  ",
                              "Note, UTF8 could cause problems with file ",
                              "names on some systems.\n");
                exit($GENERIC_ERROR);
              }
            $bchash->{$cols[0]} = join("\t",@cols[1..$num_indexes]);
          }
        else
          {print STDERR ("Unable to parse line [$ln] in barcode file ",
                         "[$bcs_file]: [",join("\t",@cols),'].  Must contain [',
                         ($num_indexes + 1),'] columns (the number of values ',
                         'supplied to --idxread plus 1 for the sample ID), ',
                         'but found [',scalar(@cols),"].  Skipping.\n")}
      }

    unless(close(FQMXBCF))
      {print STDERR ("ERROR: Unable to close fqmx barcode file: $OS_ERROR") &&
         exit($GENERIC_ERROR)}
    unless(close(BCSBCF))
      {print STDERR ("ERROR: Unable to close barcode file: $OS_ERROR") &&
         exit($OPEN_BCIN_ERROR)}

    if($dupes)
      {
        print STDERR ('ERROR: The following barcode IDs have identical index ',
                      "sequences in the barcode file: [$bcs_file]: [",
                      join(';',
                           map {join(',',@{$uniq_check->{$ARG}})}
                           grep {scalar(@{$uniq_check->{$ARG}}) > 1}
                           sort {$uniq_check->{$a}->[0] cmp
                                   $uniq_check->{$b}->[0]} keys(%$uniq_check)),
                      "].\n");
        exit($DUPE_BC_ROWS);
      }

    return($fqmx_file,$bchash,$samplehash);
  }

sub getOutfileTemplate
  {
    my $fqf_index = $ARG[0];
    my $read_num  = $fqf_index + 1;
    my $template  = '';

    if(exists($index_hash->{$fqf_index}) && !$split_all)
      {$template .= 'n/a'}
    else
      {
        if(defined($prefix))
          {$template .= $prefix}
        $template .= "\%-read-$read_num";
        if(defined($suffix))
          {$template .= $suffix}
      }

    return($template);
  }

sub getGlobPattern
  {
    my $sample       = $ARG[0];
    my $glob_pattern = '';
    if(defined($prefix))
      {$glob_pattern .= $prefix}
    $glob_pattern .= "$sample-read-*";
    if(defined($suffix))
      {$glob_pattern .= $suffix}
    return($glob_pattern);
  }

sub getDeflineSep
  {
    my $fastq_files = $ARG[0];
    my $gzipin      = $ARG[1];
    my $sep         = '';

    my $ilpat    = '[a-zA-Z0-9_:-]+ [0-9:YNATGC]+';
    my $oldilpat = '/[123fr]$';

    if($gzipin)
      {
        use IO::Uncompress::Gunzip;
        my $z = IO::Uncompress::Gunzip->new($fastq_files->[0]);
        unless($z)
          {
            print STDERR ("gunzip open failed.\n");
            exit($GENERIC_ERROR);
          }
        my $l = 0;
        while($z->getline())
          {
            $l++;
            if(/^\@/)
              {
                if(/$ilpat/)
                  {$sep = ' '}
                elsif(/$oldilpat/)
                  {$sep = '/'}
                last;
              }
            #Just in case...
            if($l > $FORMAT_DETECT_LINE_MAX)
              {last}
          }
        unless($z->close())
          {
            print STDERR ("gunzip close failed.\n");
            exit($GENERIC_ERROR);
          }
      }
    else
      {
        unless(open(FQ,$fastq_files->[0]))
          {
            print STDERR ("open of $fastq_files->[0] failed\n");
            exit($FQIN_OPEN_ERROR);
          }
        my $l = 0;
        while(<FQ>)
          {
            $l++;
            if(/^\@/)
              {
                if(/$ilpat/)
                  {$sep = ' '}
                elsif(/$oldilpat/)
                  {$sep = '/'}
                last;
              }
            #Just in case...
            if($l > $FORMAT_DETECT_LINE_MAX)
              {last}
          }
        unless(close(FQ))
          {
            print STDERR ("close of $fastq_files->[0] failed.\n");
            exit($GENERIC_ERROR);
          }
      }

    return($sep);
  }

sub fqmx2bcsStdout
  {
    my $stdout      = $ARG[0];
    my $num_indexes = $ARG[1];
    my $bcs_stdout  = '';

    unless(defined($stdout))
      {return()}

    #header
    $bcs_stdout .= join("\t",('Sample',(map {"Barcode$ARG"} (1..$num_indexes)),
                              'Count','Percent','Files')) . "\n";

    my $summary_array = [];
    my $unmatched     = join("\t",map {'unmatched'} 1..$num_indexes);
    my $count_sum     = 0;
    my $output_exists = 0;
    foreach(split(/\n/,$stdout))
      {
        #Skip header, empty, & commented lines
        if(/Id\tCount\t/ || /^\s*$/ || /^\s*#/ || /^total\t\d+$/)
          {next}

        $output_exists = 1;

        chomp;
        my @d = split(/\t/,$ARG);

        if(scalar(@d) < 3)
          {
            print STDERR ('ERROR: Unable to parse fastq-multx standard ',
                          "output.\n");
            $bcs_stdout .= "#$ARG\n";
            next;
          }

        if($d[0] ne 'unmatched' && !exists($bchash->{$d[0]}))
          {
            print STDERR ("ERROR: Unknown sample name: [$d[0]].\n");
            $bcs_stdout .= "#$ARG\n";
            next;
          }

        if($d[1] =~ /\D/)
          {
            print STDERR ("ERROR: Invalid count value: [$d[1]].\n");
            $bcs_stdout .= "#$ARG\n";
            next;
          }

        $count_sum += $d[1];
        my $glob = getGlobPattern($d[0]);
        push(@$summary_array,[$samplehash->{$d[0]},
                              ($d[0] eq 'unmatched' ?
                               $unmatched : $bchash->{$d[0]}),
                              $d[1],
                              $glob]);
      }

    foreach my $row (@$summary_array)
      {
        my $percent = ($row->[2] ? $row->[2] / $count_sum * 100 : 0);
        $bcs_stdout .= "$row->[0]\t$row->[1]\t$row->[2]\t" .
          sprintf('%.2f',$percent) . "\%\t$row->[3]\n";
      }

    return($output_exists ? $bcs_stdout : '');
  }

sub processSTDERR
  {
    my $line = $ARG[0];

    #Echo errors out to STDERR
    my $filter_pats = ['^Using Barcode File: ','End used: ',
                       #These two work around an apparently innocuous gzip issue
                       #in fastq-multx
                       'gzip: stdout: Broken pipe','^\s*$'];
    my $filter_pat  = join('|',@$filter_pats);

    #Exit non-zero when fatal error is encountered (because fastq-multx doesn't)
    my $fatal_pats =
      ['Error: number of input files \(\d+\) must match number of output files',
       'No such file or directory',
       'Returning error because of i\/o error during file close'];
    my $fatal_pat = join('|',@$fatal_pats);

    unless($line =~ /$filter_pat/i)
      {print STDERR ($line)}

    if($line =~ /$fatal_pat/)
      {return(1)}

    return(0);
  }

sub indexesFirst
  {
    my $a = $ARG[0];
    my $b = $ARG[1];

    #Determine ith both files are index reads
    my $both_r_indexes = (exists($index_hash->{$a}) &&
                          exists($index_hash->{$b}));

    #Relative order of index versus read files, index order, and read order vals
    my $idx_first = (exists($index_hash->{$b}) <=> exists($index_hash->{$a}));
    my $idx_order = ($both_r_indexes ?
                     $index_hash->{$a} <=> $index_hash->{$b} : 0);
    my $fq_order  = ($a <=> $b);

    return($idx_first || $idx_order || $fq_order);
  }

END
  {
    if(!$debug)
      {
        my $del_files = [$bcfilefqmx];
        if(defined($tmp_files))
          {push(@$del_files,@$tmp_files)}
        foreach(grep {defined($ARG) && -e $ARG} @$del_files)
          {unlink($ARG)}
      }
  }
