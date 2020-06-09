#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long qw(:config no_auto_abbrev);
use IO::Pipe::Producer;
use IO::Select;

our $VERSION     = '0.1';
our $BCSVERSION  = '0.18.6';
our $FQMXVERSION = '1.4';

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
my($bcfile,$bcfilefqmx,$stodout,$prefix,$suffix,$format,$gzipin,$deflinesep);

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
                  '<>'              => sub {push(@$fastq_files,$_[0])},
                  'fastq-multx=s'   => \$fastq_multx,
                  'debug:+'         => \$debug,
                  'check-deflines'  => \$check_deflines
                 };

processOptions();

##
## Construct the command
##

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

if($debug)
  {print STDERR ("#$command\n")}

my $producer = new IO::Pipe::Producer();
my($stdout_handle,$stderr_handle) = $producer->getSystemProducer($command);

my $sel = new IO::Select;
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
            close($fh);
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

if($? || $make_error_fatal)
  {
    print STDERR ("ERROR: fastq-multx command failed",
                  ($? ? " with a non-zero exit code [$?]" .
                   ($! =~ /./ ? " and error: $!" : '.') : ':'),
                  "\n\t$command\n");
    unless($stdout =~ /\t/)
      {exit(1)}
  }

print(fqmx2bcsStdout($stdout,scalar(@$idxread)));




sub help
  {
    my $short_only = $_[0] || 0;

    my $short_msg  = << 'END_SHORT';
usage: barcode_splitter [-h] [--version] [--bcfile FILE]
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
  -h, --help            show this help message and exit (supply twice for more)
  --version             show program's version number and exit

Barcodes:
  --bcfile FILE         REQUIRED: Tab delimited file: "Sample_ID <tab>
                        Barcode_Sequence" Multiple barcode columns with
                        different barcode lengths allowed, but all barcodes in
                        each inidividual column must be the same length.
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
END_LONG

    print($short_msg);
    print($long_msg) unless($short_only);
  }

sub help2
  {
    my $msg = << 'END_HELP2';

The following extra and modified options apply specifically to the htseq2multx.pl wrapper (as opposed to barcode_splitter).

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

sub processOptions
  {
    if(!GetOptions(%$OptionHash))
      {
        print STDERR ("Unable to parse command line.\n");
        exit(1);
      }

    if($version == 1)
      {print("$BCSVERSION\n")}
    elsif($version)
      {print("barcode_splitter, version $BCSVERSION\n",
             "htseq2multx, version $VERSION\n",
             "fastq-multx, version $FQMXVERSION\n")}
    exit(0) if($version);

    if($help == 1)
      {help() && exit(0)}
    elsif($help)
      {help2() && exit(0)}

    if(scalar(@$fastq_files) == 0)
      {
        help(1);
        print STDERR ("barcode_splitter: error: the following arguments are ",
                      "required: FILE\n");
        exit(1);
      }

    if(scalar(@$idxread) == 0)
      {
        help(1);
        print STDERR ("barcode_splitter: error: Sequence files and at least ",
                      "one number indicating the indexed file(s) (--idxread) ",
                      "is required\n");
        exit(1);
      }

    if(!defined($bcfile))
      {
        help(1);
        print STDERR ('barcode_splitter: error: Must specify a barcodes file ',
                      'with "--bcfile" option',"\n");
        exit(1);
      }
    elsif(!-e $bcfile)
      {
        #barcode_splitter does not print: help(1);
        print STDERR ("ERROR: Unable to open barcode file: [Errno 2] No such ",
                      "file or directory: '$bcfile'\n");
        exit(1);
      }

    if(scalar(grep {$_ < 1 || $_ > scalar(@$fastq_files)} @$idxread))
      {
        help(1);
        print STDERR ('barcode_splitter: error: Invalid index read number ',
                      '("--idxread"), must be between 1 and 1 (the number of ',
                      "supplied sequence files)\n");
        exit(1);
      }

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
        exit(1);
      }
    elsif(!$is_fqmx)
      {
        help2();
        print('ERROR: The fastq-multx executable appears to not be fastq-',
              'multx according to the first 2 lines of its usage.  Use the ',
              "--fastq-multx option described above to supply fastq-multx.\n");
        exit(1);
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
        exit(1);
      }

    if(!defined($prefix))
      {$prefix = ''}

    if(defined($prefix) && $galaxy)
      {$prefix =~ s/_/-/g}

    if(defined($suffix) && $galaxy)
      {$suffix =~ s/_/-/g}

    #Create a hash of the index reads
    my $order = 1;
    foreach my $fqf_num (@$idxread)
      {$index_hash->{$fqf_num - 1} = $order++}

    ($bcfilefqmx,$bchash,$samplehash) =
      bcs2fqmxBCFile($bcfile,scalar(@$idxread));

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

    if(!$split_all && scalar(@$fastq_files) == scalar(@$idxread))
      {$split_all = 1}

    if(outfilesExist($bchash,$fastq_files,$split_all))
      {exit(7)}

    if(isOutdirMissing($prefix))
      {
        print STDERR ("ERROR: Directory in prefix [$prefix] does not exist.\n");
        exit(1);
      }

    my $missing_infiles = [grep {!-e $_} @$fastq_files];
    if(scalar(@$missing_infiles))
      {
        print STDERR ('ERROR: ',scalar(@$missing_infiles),' input files do ',
                      'not exist: ',join(',',@$missing_infiles),"\n");
        exit(1);
      }
  }

#This detects and converts files input via process substitution, and writes them
#out to temporary files.  Assumes $gzip is true.
sub processGZInfiles
  {
    my $fastq_files = $_[0];
    my $prefix      = $_[1];
    return() if(scalar(grep {$_ !~ /\.gz$/i} @$fastq_files) == 0);

    my $tmp_files = [];

    #Create a temporary set of files
    my $time = time();
    my $dir = '.';
    if(defined($prefix) && $prefix =~ m%(.*/)%)
      {$dir = $1}
    my $template = "$dir/${time}-\%.fq.gz";

    for(my $i = 0;$i <= $#{$fastq_files};$i++)
      {
        next if($fastq_files->[$i] !~ /\.gz$/i);
        my $file = $template;
        $file =~ s/%/$i/;
        $fastq_files->[$i] = $file;
        push(@$tmp_files,$file);
        writeGZProcSubFile($fastq_files->[$i],$file);
      }

    return(wantarray ? @$tmp_files : $tmp_files);
  }

sub writeGZProcSubFile
  {
    my $inprocsub = $_[0];
    my $outfile   = $_[1];

    unless(open(IN,$inprocsub))
      {
        print STDERR ('ERROR: Unable to open input stream from process ',
                      "substitution [$inprocsub]\n");
        exit(1);
      }
    binmode(IN);
    unless(open(OUT,'>',$outfile))
      {
        print STDERR ('ERROR: Unable to open temporary output file [$outfile]',
                      " for process substitution [$inprocsub].\n");
        exit(1);
      }
    binmode(OUT);

    print OUT (<IN>);

    close(OUT);
    close(IN);
  }

sub isOutdirMissing
  {
    my $prefix  = $_[0];
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
    my $hash     = $_[0];
    my $files    = $_[1];
    my $existing = 0;

    foreach my $fqf_index (0..$#{$files})
      {
        my $template = getOutfileTemplate($fqf_index);
        next if($template =~ m%n/a%i);
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

    my $output  = `($fastq_multx 3>&1 1>&2- 2>&3- ) | head -n 2`;
    if($?)
      {
        $error = "fastq-multx failed with a non-zero exit code [$?]" .
          ($! =~ /./ ? " and error: $!" : '.') . "\n";
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
    my $used = $_[0];
    my $reqd = $_[1];

    my $lim =
      (scalar(@$used) < scalar(@$reqd) ? scalar(@$used) : scalar(@$reqd)) - 1;

    return(1) if($lim == 0);

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

sub bcs2fqmxBCFile
  {
    my $bcs_file    = $_[0];
    my $num_indexes = $_[1];
    my $fqmx_file   = $bcs_file . '.fqmx';
    my $bchash      = {};
    my $samplehash  = {unmatched => 'unmatched'};
    my $uniq_check  = {};
    my $dupes       = 0;

    unless(open(BCSBCF,$bcs_file))
      {print STDERR ("ERROR: Unable to open barcode file: $!") && exit(1)}

    unless(open(FQMXBCF,'>',$fqmx_file))
      {print STDERR ("ERROR: Unable to open fqmx barcode file: $!") && exit(1)}

    my $ln     = 0;
    while(<BCSBCF>)
      {
        $ln++;
        chomp;
        next if(/^\s*#/ || /^\s*$/);
        my @cols = split(/\t/,$_);
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
                exit(1);
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

    close(FQMXBCF);
    close(BCSBCF);

    if($dupes)
      {
        print STDERR ('ERROR: The following barcode IDs have identical index ',
                      "sequences in the barcode file: [$bcs_file]: [",
                      join(';',
                           map {join(',',@{$uniq_check->{$_}})}
                           grep {scalar(@{$uniq_check->{$_}}) > 1}
                           sort {$uniq_check->{$a}->[0] cmp
                                   $uniq_check->{$b}->[0]} keys(%$uniq_check)),
                      "].\n");
        exit(12);
      }

    return($fqmx_file,$bchash,$samplehash);
  }

sub getOutfileTemplate
  {
    my $fqf_index = $_[0];
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
    my $sample       = $_[0];
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
    my $fastq_files = $_[0];
    my $gzipin      = $_[1];
    my $sep         = '';

    my $ilpat    = '[a-zA-Z0-9_:-]+ [0-9:YNATGC]+';
    my $oldilpat = '/[123fr]$';

    if($gzipin)
      {
        use IO::Uncompress::Gunzip;
        my $z = new IO::Uncompress::Gunzip $fastq_files->[0] ||
          die("gunzip failed\n");
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
            last if($l > 100);
          }
        $z->close();
      }
    else
      {
        unless(open(FQ,$fastq_files->[0]))
          {die("open failed\n")}
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
            last if($l > 100);
          }
      }

    return($sep);
  }

sub fqmx2bcsStdout
  {
    my $stdout      = $_[0];
    my $num_indexes = $_[1];
    my $bcs_stdout  = '';

    return() unless(defined($stdout));

    #header
    $bcs_stdout .= join("\t",('Sample',(map {"Barcode$_"} (1..$num_indexes)),
                              'Count','Percent','Files')) . "\n";

    my $summary_array = [];
    my $unmatched     = join("\t",map {'unmatched'} 1..$num_indexes);
    my $count_sum     = 0;
    my $output_exists = 0;
    foreach(split(/\n/,$stdout))
      {
        #Skip header, empty, & commented lines
        next if(/Id\tCount\t/ || /^\s*$/ || /^\s*#/ || /^total\t\d+$/);

        $output_exists = 1;

        chomp;
        my @d = split(/\t/,$_);

        if(scalar(@d) < 3)
          {
            print STDERR ('ERROR: Unable to parse fastq-multx standard ',
                          "output.\n");
            $bcs_stdout .= "#$_\n";
            next;
          }

        if($d[0] ne 'unmatched' && !exists($bchash->{$d[0]}))
          {
            print STDERR ("ERROR: Unknown sample name: [$d[0]].\n");
            $bcs_stdout .= "#$_\n";
            next;
          }

        if($d[1] =~ /\D/)
          {
            print STDERR ("ERROR: Invalid count value: [$d[1]].\n");
            $bcs_stdout .= "#$_\n";
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
    my $line = $_[0];

    #Echo errors out to STDERR
    my $filter_pats = ['^Using Barcode File: ','End used: '];
    my $filter_pat  = join('|',@$filter_pats);

    #Exit non-zero when fatal error is encountered (because fastq-multx doesn't)
    my $fatal_pats  = ['Error: number of input files \(\d+\) must match ' .
                       'number of output files',
                       'No such file or directory'];
    my $fatal_pat   = join('|',@$fatal_pats);

    unless($line =~ /$filter_pat/i)
      {print STDERR ($line)}

    if($line =~ /$fatal_pat/)
      {return(1)}

    return(0);
  }

sub indexesFirst
  {
    my $a = $_[0];
    my $b = $_[1];

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
        push(@$del_files,@$tmp_files) if(defined($tmp_files));
        foreach(grep {defined($_) && -e $_} @$del_files)
          {unlink($_)}
      }
  }
