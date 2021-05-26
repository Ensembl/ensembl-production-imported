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

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $registry = 'Bio::EnsEMBL::Registry';

my $reg_path = $opt{registry};
$registry->load_all($reg_path);

my $sps = $registry->get_all_species();

say scalar(@$sps) . " species";
for my $sp (sort @$sps) {

  my $dbas = $registry->get_all_DBAdaptors($sp);
  my %groups = map { $_->group => 1 } @$dbas;
  
  my $stats = "";
  my $db = "";
  my $name = "";
  my ($core) = grep { $_->group eq 'core' } @$dbas;
  if ($core) {
    try {
      $db = $core->dbc->dbname;
      my $genea = $core->get_GeneAdaptor();
      my $tra = $core->get_TranscriptAdaptor();
      my $meta = $registry->get_adaptor($sp, "core", "MetaContainer");
      my ($insdc) = @{ $meta->list_value_by_key("assembly.accession") };
      $stats .= "$insdc\t" if $insdc;
      
      # BRC4 specific
      my ($org) = @{ $meta->list_value_by_key("BRC4.organism_abbrev") };
      my ($comp) = @{ $meta->list_value_by_key("BRC4.component") };
      $name = "$org\t" if $org;
      $name .= "$comp\t" if $comp;
      
      $core->dbc->disconnect_if_idle();
    } catch {
      warn("Error: can't use core for $sp");
    };
  }

  say "$db\t$sp\t$name\t" . join(", ", sort keys %groups) . "\t$stats";
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
    Show species in a registry

    --registry <path> : Ensembl registry
    
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
    "help",
    "verbose",
    "debug",
  );

  usage("Registry needed") if not $opt{registry};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

