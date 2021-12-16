#!/usr/env perl
use v5.14.00;
use strict;
use warnings;
use Carp;
use autodie qw(:all);
use Readonly;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl qw( :easy ); 
Log::Log4perl->easy_init($WARN); 
my $logger = get_logger(); 

use File::Spec::Functions qw(catdir catfile);
use JSON;
use Data::Dumper;
use List::MoreUtils qw(zip);
use Try::Tiny;

###############################################################################
# MAIN
# Get command line args
our %opt = %{ opt_check() };

my $input = get_input_metadata($opt{input});
get_dir_data($opt{dir}, $input);

sub get_dir_data {
  my ($dir, $input) = @_;

  my @header = qw(Component Species Dataset Sample Mapped Reads InputPaired InputStranded Paired Stranded DiffPaired DiffStranded);
  say join("\t", @header);
  for my $component (sub_dirs($dir)) {
    my $comp_dir = catdir($dir, $component);
    for my $species (sub_dirs($comp_dir)) {
      my $species_dir = catdir($comp_dir, $species);
      for my $dataset (sub_dirs($species_dir)) {
        my $dataset_dir = catdir($species_dir, $dataset);
        my $dataset_stranded = 0;
        my $dataset_unstranded = 0;
        for my $sample (sub_dirs($dataset_dir)) {
          my %data = (
            Component => $component,
            Species => $species,
            Dataset => $dataset,
            Sample => $sample,
          );
          my $sample_dir = catdir($dataset_dir, $sample);
          
          # Get sample metadata
          my $sample_file = catfile($sample_dir, 'metadata.json');
          my $metadata = get_sample_metadata($sample_file);
          $data{Paired} = $metadata->{hasPairedEnds} ? 1 : 0;
          $data{Stranded} = $metadata->{isStrandSpecific} ? 1 : 0;
          $dataset_stranded++ if $data{Stranded};
          $dataset_unstranded++ if not $data{Stranded};

          # Get results stats
          my $results_file = catfile($sample_dir, 'mappingStats.txt');
          my $res_stats = get_results_stats($results_file);
          $data{Mapped}= $res_stats->{mapped} ? $res_stats->{mapped} : "NA";
          $data{Reads}= $res_stats->{number_reads_mapped} ? $res_stats->{number_reads_mapped} : "NA";
          
          my $input_metadata = $input->{$component}->{$species}->{$dataset}->{$sample};
          if ($input_metadata) {
            $data{InputPaired} = $input_metadata->{hasPairedEnds} ? 1 : 0;
            $data{InputStranded} = $input_metadata->{isStrandSpecific} ? 1 : 0;
            
            if ($data{Paired} != $data{InputPaired}) {
              say STDERR "PAIRED DIFF\t$component\t$species\t$dataset\t$sample" if $opt{v};
              $data{DiffPaired} = 1;
            } else {
              $data{DiffPaired} = 0;
            }
            if ($data{Stranded} != $data{InputStranded}) {
              say STDERR "STRAND DIFF\t$component\t$species\t$dataset\t$sample" if $opt{v};
              $data{DiffStranded} = 1;
            } else {
              $data{DiffStranded} = 0;
            }
          } else {
            say STDERR "{$component}->{$species}->{$dataset}->{$sample} has no metadata?";
            next;
          }
          
          say join("\t", map { exists $data{$_} ? $data{$_} : "NA" } @header);
        }
        
        # Check the dataset is homogenous
        if ($dataset_stranded and $dataset_unstranded) {
          say STDERR "STRAND MIX\t$component\t$species\t$dataset";
        }
      }
    }
  }
}

sub sub_dirs {
  my ($dir) = @_;
  
  opendir(my $dh, $dir);
  my @sub_dirs = grep { not /^\./ and -d catdir($dir, $_) } readdir $dh;
  closedir $dh; 
  return @sub_dirs;
}

sub get_sample_metadata {
  my ($path) = @_;
  
  if (-s $path) {
    return read_json($path);
  } else {
    return {};
  }
}

sub read_json {
  my ($path) = @_;
  
  my $json_str;
  open my $fh, "<", $path;
  while (my $line = readline $fh) {
    chomp $line;
    $json_str .= $line;
  }
  close $fh;
  
  return decode_json($json_str);
}

sub get_results_stats {
  my ($path) = @_;
  
  my %stats;
  try {
    open my $infh, "<", $path;
    my $first_line = readline $infh;
    chomp $first_line;
    my @header = split "\t", $first_line;

    while (my $line = readline $infh) {
      chomp $line;
      if ($line =~ /^results.bam/) {
        my @data = split "\t", $line;
        %stats = zip @header, @data;
        last;
      }
    }
    close $infh;
  } catch {
    say STDERR "Can't read $path";
  };
  return \%stats;
}

sub get_input_metadata {
  my ($path) = @_;
  
  my $datasets = read_json($path);
  
  my %input;
  for my $dataset (@$datasets) {
    my $component = $dataset->{component};
    my $species = $dataset->{species};
    my $dataset_name = $dataset->{name};
    for my $sample (@{ $dataset->{runs} }) {
      my $name = $sample->{name} // $sample->{accessions}->[0];
      $input{$component}{$species}{$dataset_name}{$name} = $sample;
    }
  }
  
  return \%input;
}

###############################################################################
# Parameters and usage
sub usage {
  my $error = shift;
  my $help = '';
  if ($error) {
    $help = "[ $error ]\n";
  }
  $help .= <<'EOF';
    Write a summary of files generated by the BRC4 RNA-Seq pipeline
    
    --dir <path>      : Path to the root of the RNA-Seq data
    --input <path>    : Input json file from Eupath with all datasets
    
    --help            : show this help message
    --verbose         : show detailed progress
    --debug           : show even more information (for debugging purposes)
EOF
  print STDERR "$help\n";
  exit(1);
}

sub opt_check {
  my %opt = ();
  GetOptions(\%opt,
    "dir=s",
    "input=s",
    "help",
    "verbose",
    "debug",
  );

  usage()                if $opt{help};
  usage("Dir needed")    if not $opt{dir};
  usage("Input json needed")    if not $opt{input};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}
__END__


