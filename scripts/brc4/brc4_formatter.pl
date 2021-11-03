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


my $species = $opt{species};

my $component = $opt{component};
my $organism = $opt{organism};
my $update = $opt{update};

format_brc4_db($registry, $species, $component, $organism, $update);

###############################################################################

sub format_brc4_db {
  my ($registry, $species, $component, $organism, $update) = @_;
  
  say STDERR "Format database for $species to BRC4 specs";

  my $meta = $registry->get_adaptor($species, "core", "MetaContainer");
  my $sa = $registry->get_adaptor($species, "core", "Slice");
  my $ata = $registry->get_adaptor($species, 'core', 'attribute');
  
  # Update component and organism abbrev
  add_meta_value($meta, 'BRC4.component', $component, $update);
  add_meta_value($meta, 'BRC4.organism_abbrev', $organism, $update);
  
  # Add seq_region_attribs
  add_seqregions_meta($sa, $ata, $update);
  
  # Add gene, transcript, translation, exon versions if there are none
  add_versions($sa->dbc, $update);
}

sub add_versions {
  my ($dbc, $update) = @_;
  
  for my $table (qw{ gene transcript translation exon }) {
    my $count_query = "SELECT * FROM $table WHERE version IS NULL";
    my $sth = $dbc->prepare($count_query);
    $sth->execute();
    my $count = $sth->rows;
    say STDERR "$count features without version in table $table";
    if ($update) {
      my $update_query = "UPDATE $table SET version = 1 WHERE version IS NULL";
      $dbc->do($update_query);
    }
  }
}

sub add_meta_value {
  my ($meta, $key, $new_value, $update) = @_;
  
  my ($cur_value) = @{ $meta->list_value_by_key($key) };
  
  if ($cur_value) {
    say STDERR "WARNING: there is already a value for meta key '$key': '$cur_value'";
  } else {
    if ($update) {
      say STDERR "Insert new value '$new_value' for key '$key'";
      # TODO: actually insert the meta data
      my $add_query = "INSERT INTO meta(meta_key, meta_value) VALUES('$key', '$new_value')";
      $meta->dbc->do($add_query);
    } else {
      say STDERR "Dry run: no insertion of new_value for $key";
    }
  }
}

sub add_seqregions_meta {
  my ($sa, $ata, $update) = @_;
  
  my $slices = $sa->fetch_all('toplevel');
  
  for my $slice (@$slices) {
    my $seq_name = $slice->seq_region_name;
    
    # Check their attribs and synonyms
    my $attribs = $slice->get_all_Attributes();
    my $synonyms = $slice->get_all_synonyms();
    
    my ($cur_ebi_name, $cur_brc4_name);
    for my $at (@$attribs) {
      $cur_ebi_name = $at->value if $at->code eq 'EBI_seq_region_name';
      $cur_brc4_name = $at->value if $at->code eq 'BRC4_seq_region_name';
    }
    
    my ($refseq_name, $insdc_name);
    for my $syn (@$synonyms) {
      $insdc_name = $syn->name if $syn->dbname eq 'INSDC';
      $refseq_name = $syn->name if $syn->dbname eq 'RefSeq';
    }
    my $new_brc4_name = $insdc_name // $refseq_name;
    
    # Check EBI name
    if ($cur_ebi_name) {
      if ($cur_ebi_name ne $seq_name) {
        say STDERR "WARNING: current EBI_seq_region_name does not match the seq_region name: $cur_ebi_name vs $seq_name";
      } else {
        say "The EBI name for $seq_name is already properly defined";
      }
    } else {
      if ($update) {
        say STDERR "Insert new EBI_seq_region_name as $seq_name";
        insert_attrib($ata, $slice, 'EBI_seq_region_name', $seq_name);
      }
    }

    # Check BRC4 name
    if ($cur_brc4_name) {
      if ($cur_brc4_name ne $new_brc4_name) {
        say STDERR "WARNING: current BRC_seq_region_name does not match the expected seq_region name: $cur_brc4_name vs $new_brc4_name";
      } else {
        say "The BRC4 name for $seq_name is already properly defined as $new_brc4_name";
      }
    } else {
      if ($update) {
        say STDERR "Insert new BRC4_seq_region_name as $new_brc4_name";
        insert_attrib($ata, $slice, 'BRC4_seq_region_name', $new_brc4_name);
      }
    }
  }
}

sub insert_attrib {
  my ($ata, $slice, $name, $value) = @_;

  my $attr = Bio::EnsEMBL::Attribute->new(
    -CODE => $name,
    -VALUE => $value,
  );
  my @attrs = ($attr);

  $ata->store_on_Slice($slice, \@attrs);
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
    Format a database from Ensembl to follow BRC4 specifications.

    --registry <path> : Ensembl registry
    --species <str>   : production_name of one species

    --component <str> : BRC4 component for that species
    --organism <str>  : BRC4 organism_abbrev for that species
    
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
    "component=s",
    "organism=s",
    "update",
    "help",
    "verbose",
    "debug",
  ) or usage();

  usage("Registry needed") if not $opt{registry};
  usage("Species needed") if not $opt{species};
  usage()                if $opt{help};
  
  if ($opt{update} and not ($opt{component} and $opt{organism})) {
    usage("Can't update the database without a component and organism");
  }
  
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

