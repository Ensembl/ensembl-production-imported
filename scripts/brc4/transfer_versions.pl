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

my %old_data;
for my $feat (@features) {
  $old_data{$feat} = get_old_data($registry, $species, $feat);
  print STDERR scalar(keys $old_data{$feat}) . " ${feat}s from old database\n";
}

# Reload registry, this time for new database
$reg_path = $opt{new_registry};
$registry->load_all($reg_path);

if ($opt{descriptions}) {
  say STDERR "Gene descriptions transfer:";
  update_descriptions($registry, $species, $old_data{gene}, $opt{write});
}

if ($opt{versions}) {
  say STDERR "Gene versions transfer:";
  # Update the version for the features we want to transfer and update
  for my $feat (@features) {
    update_versions($registry, $species, $feat, $old_data{$feat}, $opt{write});
  }

  # Blanket version replacement for transcripts, translations, exons
  if ($opt{write}) {
    say STDERR "Transcripts, translations and exons version update";
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
    say STDERR "Transcripts, translations and exons init versions were not updated (use --write to do so)";
  }
}

###############################################################################
sub get_old_data {
  my ($registry, $species, $feature) = @_;
  
  my %feats;
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version;
    my $description;
    if ($feature eq 'gene') {
      $description = $feat->description;
    }
    $feats{$id} = { version => $version, description => $description };
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
      my $old_feat = $old_feats->{$id};
      
      if (not $old_feat) {
        $new_count++;
        next;
      }
      
      my $old_version = $old_feat->{version};
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
  print STDERR "$new_count $feature version to be initialized\n";
  print STDERR "(Use --write to make the changes to the database)\n" if $update_count + $new_count > 0 and not $update;
}

sub update_descriptions {
  my ($registry, $species, $old_genes, $update) = @_;
  
  my $update_count = 0;
  my $empty_count = 0;
  my $new_count = 0;
  
  my $ga = $registry->get_adaptor($species, "core", 'gene');
  for my $gene (@{$ga->fetch_all}) {
    my $id = $gene->stable_id;
    my $description = $gene->description;

    if (not defined $description) {
      my $old_gene = $old_genes->{$id};
      if (not $old_gene) {
        $new_count++;
        next;
      }
      my $old_description = $old_gene->{description};

      if (defined $old_description) {
        my $new_description = $old_description;
        say STDERR "Transfer gene $id description: $new_description" if $opt{verbose};
        $update_count++;

        if ($update) {
          $gene->description($new_description);
          $ga->update($gene);
        }
      } else {
        $empty_count++;
      }
    }
  }
  
  print STDERR "$update_count gene descriptions transferred\n";
  print STDERR "$empty_count genes without description remain\n";
  print STDERR "$empty_count new genes, without description\n";
  print STDERR "(Use --write to update the descriptions in the database)\n" if $update_count > 0 and not $update;
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
    Transfer the gene versions and descriptions to a patched database

    --old_registry <path> : Ensembl registry
    --new_registry <path> : Ensembl registry
    --species <str>   : production_name of one species

    --descriptions    : Transfer the gene descriptions
    --versions        : Transfer the gene versions, and init the others
    
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
    "descriptions",
    "versions",
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

