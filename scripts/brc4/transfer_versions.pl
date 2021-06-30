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
our %opt = %{ opt_check() };

my $species = $opt{species};

# Features to transfer versions from
my @features = qw(gene);

# Get genes and transcripts from old database
my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{old_registry};
$registry->load_all($reg_path);

my %old_versions;
for my $feat (@features) {
  $old_versions{$feat} = get_versions($registry, $species, $feat);
  print STDERR scalar(keys $old_versions{$feat}) . " ${feat}s from old database\n";
}

# Reload registry, this time for new database
$reg_path = $opt{new_registry};
$registry->load_all($reg_path);

# Update the version for the features we want to transfer and update
for my $feat (@features) {
  update_versions($registry, $species, $feat, $old_versions{$feat}, $opt{write});
}

# Blanket version replacement for transcripts, translations, exons
if ($opt{write}) {
  my $ga = $registry->get_adaptor($species, "core", 'gene');
  my $dbc = $ga->dbc;
  
  # Transfer the versions from genes to transcripts
  my $transcript_query = "UPDATE transcript LEFT JOIN gene USING(gene_id) SET transcript.version = gene.version";
  $dbc->do($transcript_query);

  # Also transfer the versions from transcripts to translations
  my $translation_query = "UPDATE translation LEFT JOIN transcript USING(transcript_id) SET translation.version = transcript.version";
  $dbc->do($translation_query);

  # And set the exons version to 1
  my $exon_query = "UPDATE exon SET version = 1";
  $dbc->do($exon_query);
} else {
  say STDERR "Transcripts, translations and exons were not updated (use --write to do so)";
}

###############################################################################
sub get_versions {
  my ($registry, $species, $feature) = @_;
  
  my %feats;
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version;
    $feats{$id} = $version;
  }
  
  return \%feats;
}

sub update_versions {
  my ($registry, $species, $feature, $old_feats, $update) = @_;
  
  my $update_count = 0;
  my $new_count = 0;
  
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version;

    if (not defined $version) {
      my $old_version = $old_feats->{$id};
      my $new_version = 1;

      if (defined $old_version) {
        $new_version = $old_version + 1;
        say STDERR "Updated $feature $id must be upgrade from $old_version to $new_version" if $opt{verbose};
        $update_count++;
      } else {
        say STDERR "New $feature $id must be initialized to 1" if $opt{verbose};
        $new_count++;
      }

      if ($update) {
        $feat->version($new_version);
        $fa->update($feat);
      }
    }
  }
  
  print STDERR "$update_count $feature version updated\n";
  print STDERR "$new_count $feature version initialized\n";
  print STDERR "(Use --write to make the changes to the database)\n" if $update_count + $new_count > 0 and not $update;
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
    Transfer and update genes and transcripts versions to a patched database

    --old_registry <path> : Ensembl registry
    --new_registry <path> : Ensembl registry
    --species <str>   : production_name of one species
    
    --write           : Do the actual changes (default is no changes to the database)
    
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
    "old_registry=s",
    "new_registry=s",
    "species=s",
    "write",
    "help",
    "verbose",
    "debug",
  );

  usage("Old registry needed") if not $opt{old_registry};
  usage("New registry needed") if not $opt{new_registry};
  usage("Species needed") if not $opt{species};

  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

