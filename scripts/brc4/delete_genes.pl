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

use Bio::EnsEMBL::Registry;
use Try::Tiny;
use File::Path qw(make_path);

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{registry};
$registry->load_all($reg_path);

my $list = get_genes($opt{list});
print STDERR scalar(@$list) . " genes to delete\n";
delete_genes($registry, $opt{species}, $list, $opt{update});

sub get_genes {
  my ($path) = @_;
  
  my @genes;
  
  open my $fh, "<", $path;
  while(my $line = readline $fh) {
    chomp $line;
    $line =~ s/\s+//g;
    if ($line) {
      push @genes, $line;
    }
  }
  close $fh;
  
  return \@genes;
}

sub delete_genes {
  my ($registry, $species, $list, $update) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  
  my $delete_count = 0;
  
  for my $id (@$list) {
    my $gene = $ga->fetch_by_stable_id($id);
    
    if (not $gene) {
      print STDERR "WARNING: No gene found with stable_id = '$id'\n";
      next;
    }
    
    if ($update) {
      $ga->remove($gene);
      $delete_count++;
    }
  }
  
  print STDERR "$delete_count genes deleted\n";
  print STDERR "(Use --update to make the changes to the database)\n" if $delete_count == 0 and @$list > 0;
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
    Delete a list of genes from a core db

    --registry <path> : Ensembl registry
    --species <str>   : production_name of one species
    --list <path>     : List of genes to delete
    
    --update          : Do the actual deletion (default is no changes to the database)
    
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
    "registry=s",
    "species=s",
    "list=s",
    "update",
    "help",
    "verbose",
    "debug",
  );

  usage("Registry needed") if not $opt{registry};
  usage("Species needed") if not $opt{species};
  usage("List needed") if not $opt{list};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

