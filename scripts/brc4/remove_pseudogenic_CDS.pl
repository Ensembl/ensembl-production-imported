#!/usr/env perl
=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


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

use Data::Dumper;

###############################################################################
# MAIN
main();

sub main {
  # Get command line args
  my %opt = %{ opt_check() };

  $logger->info("Load registry");
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($opt{registry}, 1);

  my @all_species = ($opt{species}) || @{$registry->get_all_species()};
  for my $species (sort @all_species) {
    check_pseudoCDS($registry, $species, $opt{update});
  }
}

sub check_pseudoCDS {
  my ($registry, $species, $update) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "Gene");
  my $tra = $registry->get_adaptor($species, "core", "Transcript");
  my $tla = $registry->get_adaptor($species, "core", "Translation");
  my $dbname = $tra->dbc->dbname;
  $logger->info("Database:\t$dbname");
  $logger->info("Species:\t$species");
  
  # We're looking at all pseudogenes/+with_CDS
  # as well as protein_coding genes, in case they contain both protein_coding and
  # pseudogenic transcripts
  my $all_genes = get_genes($tra);
  my $count_pst = 0;

  GENE: for my $gene (@$all_genes) {
    my $not_pseudo = 0;

    for my $tr (@{$gene->get_all_Transcripts()}) {
      my $translation = $tr->translation;
      my $tr_biotype = $tr->biotype;

      # Only consider transcripts that are exactly pseudogene_with_CDS,
      # even if the gene is something else
      if ($tr_biotype eq 'pseudogene_with_CDS') {
        $count_pst++;
        
        # Update the transcript biotype to pseudogene + remove its translation
        if ($update) {
          $logger->debug("UPDATE transcript " . $tr->stable_id);
          $tla->remove($translation);
          $tr->biotype('pseudogene');
          $tra->update($tr);
        }
      # Note if the gene contains transcripts that are not pseudogenic
      } else {
        if ($tr_biotype ne 'pseudogene') {
          $not_pseudo++;
        }
      }
    }

    # Update the gene biotype to pseudogene too, but only if all transcripts are pseudogenic
    my $all_pseudo = not $not_pseudo;
    if ($update and $all_pseudo) {
      $logger->debug("UPDATE gene " . $gene->stable_id);
      $gene->biotype('pseudogene');
      $ga->update($gene);
    }
  }
  $ga->dbc->disconnect_if_idle();
  $logger->info("Genes with pseudogene_with_CDS transcripts\t" . scalar(@$all_genes));
  $logger->info("Pseudogenic transcript with CDS\t$count_pst");
}

sub get_genes {
  my ($tra) = @_;

  # Get the genes that contain at least one transcript with biotype pseudogene_with_CDS
  my @genes;
  for my $tr (@{$tra->fetch_all_by_biotype('pseudogene_with_CDS')}) {
    push @genes, $tr->get_Gene();
  }
  return \@genes;
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
    Check and remove pseudogenic_CDSs

    --registry <path> : Ensembl registry for the core database

    Optional:
    --species <str>   : production_name of one species
    --update          : remove the pseudogenic_CDSs
    
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
    "update",
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

