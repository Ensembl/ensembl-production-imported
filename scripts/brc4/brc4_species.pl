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
use Capture::Tiny ':all';

use Bio::EnsEMBL::Registry;
use Try::Tiny;

# Ordered list of meta fields to use
my @fields = qw(
BRC4.component
species.scientific_name
species.strain
BRC4.organism_abbrev
assembly.accession
);

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $registry = 'Bio::EnsEMBL::Registry';

my $reg_path = $opt{registry};
$registry->load_all($reg_path);

my $sps = $registry->get_all_species();

# Hash to track unique key-value pairs
my %unique_abbrevs;

my @genomes;
for my $sp (sort @$sps) {
  my $dbas;
  my %groups;
  $dbas = $registry->get_all_DBAdaptors($sp);
  %groups = map { $_->group => 1 } @$dbas;
  
  my $db = "";
  my $name = "";
  my ($core) = grep { $_->group eq 'core' } @$dbas;
  my $skip = 0;
  my %stats;

  if ($core) {
    try {
      my ($stdout, $stderr) = capture {
        $db = $core->dbc->dbname;
        my $genea = $core->get_GeneAdaptor();
        my $tra = $core->get_TranscriptAdaptor();
        my $meta = $registry->get_adaptor($sp, "core", "MetaContainer");
        
        for my $key (@fields) {
          $stats{$key} = get_meta_value($meta, $key);
        }
      };
      $core->dbc->disconnect_if_idle();
      print($stdout);
      
      print STDERR $stderr if $opt{debug};
    } catch {
      warn("Error: can't use core for $sp: $_");
    };
  }

  # To print
  push @genomes, \%stats;
}

# Print all genomes metadata in order
for my $genome (sort {
    $a->{'BRC4.component'} cmp $b->{'BRC4.component'}
      or $a->{'species.scientific_name'} cmp $b->{'species.scientific_name'}
      or $a->{'BRC4.organism_abbrev'} cmp $b->{'BRC4.organism_abbrev'}
  } @genomes) {
  # Check if BRC4.organism_abbrev is unique
  my $abbrev = $genome->{'BRC4.organism_abbrev'};
  my $prod_name = $genome->{'species.scientific_name'};

  if (defined $abbrev && !$unique_abbrevs{$abbrev}) {
      # BRC4.organism_abbrev is unique
      $unique_abbrevs{$abbrev}=$prod_name;
  }
  else{
    die "Error: Non-unique abbreviation encountered: $abbrev $prod_name\n";
  }
  say join("\t", map { $genome->{$_} // "" } @fields);
  }

# Return the value for all the keys
sub get_meta_value {
  my ($meta, $key) = @_;
  my ($value) =@{ $meta->list_value_by_key($key) };
  return $value;
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
    Create a list of all species metadata in BRC4 prod.

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
