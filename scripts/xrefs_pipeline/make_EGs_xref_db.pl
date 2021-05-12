#!/usr/bin/env perl -w
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# All EG specific Display Xrefs, like JGI or Broad or GeneDB ones

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

my $regfile = $ENV{'REGISTRY'};
my $user = $ENV{'USER'};
my $species = '';
my $all = '';
my $help;
my $out;
$| = 1;

GetOptions(
     'regfile:s'  => \$regfile,
     'species:s'  => \$species,
     'all!'       => \$all,
     'help!'      => \$help,
     'h!'         => \$help,     
     'out=s'      => \$out,
);

die("Cannot do anything without a kill list file and a registry file and a compara database
  'regfile:s'   => $regfile,
  'species:s'  => $species, (single species to run on )
  'all!'       => $all, ( run on all species )
") if ($species eq '' and $all eq '' or $help);

Bio::EnsEMBL::Registry->load_all($regfile);

open ( SQL,">$out") or die("Cannot open file $out for writing\n");

my @dbas;

if ($all){
  @dbas = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors()};
}

if ($species ne ''){
  my @db = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors(-species => "$species")};
  die ("Cannot find adaptor for species $species in registry\n") unless scalar(@db) == 1;
  push @dbas, @db;
}

@dbas = sort { $a->species cmp $b->species } @dbas;

foreach my $dba (@dbas){

    # print STDERR "db name: " . $dba->dbc->dbname . "\n";
    
    next unless $dba->dbc->dbname !~ /collection/;
    
    next unless $dba->group eq 'core';
    if ( $dba->species =~ /(\w+)\s(\w+)/ ){
  my $shortname = uc(substr($1,0,4));
  $shortname .= "_". uc(substr($2,0,1));
  my $name = $dba->species;
  print "NAME $name\n";
  my $meta = $dba->get_MetaContainer;
  my $tax_id = $meta->get_taxonomy_id;
  
  # we want to fetch all the protein coding genes and their associated xrefs
  my $ga = $dba->get_GeneAdaptor;
  
  print STDERR "got " . @{$ga->fetch_all_by_logic_name('protein_coding')} . " protein coding genes\n";
  
  my @protein_coding_genes = @{$ga->fetch_all_by_biotype('protein_coding')};
  push (@protein_coding_genes, @{$ga->fetch_all_by_biotype('pseudogene')});
  push (@protein_coding_genes, @{$ga->fetch_all_by_biotype('non_translating_cds')});
  
  foreach my $pc_gene ( @protein_coding_genes ) {
      
      foreach my $pcRNA ( @{$pc_gene->get_all_Transcripts} ) {
    
    # Get the display_xref and keep it if it is ENA_GENE or PGD_GENE

    my $display_db_entry_xref = $pc_gene->display_xref;
    if ($display_db_entry_xref->dbname =~ /ena_gene|pgd_gene/i) {
        print SQL "INSERT INTO EG_Xref values (\\N,\"" . 
      $name . '",' .
      $tax_id. ',"' .
      $pc_gene->stable_id . '",' .
      $pc_gene->dbID . ',"' .
      $pcRNA->stable_id .'",' .
      $pcRNA->dbID .',"' .
      $pcRNA->biotype . '","' .
      $display_db_entry_xref->database .'","' .
      $display_db_entry_xref->display_id .'","' .
      $display_db_entry_xref->primary_id .'","' .
      $display_db_entry_xref->description . '");'."\n" ;
    }
    
    # Xrefs are at the Gene levels in our case !
    
    foreach my $xref ( @{$pc_gene->get_all_DBEntries} ) {
        next unless $xref->database =~ /broad/i or $xref->database =~ /jgi/i 
          or $xref->database =~ /genedb/i or $xref->database =~ /cadre/i 
          or $xref->database =~ /aspgd/i
          or $xref->database =~ /ena_gene/i or $xref->database =~ /pgd_gene/i
          or $xref->database =~ /schistodb/i or $xref->database =~ /mycgr3_jgi_v2.0_gene/i;
        my $description = $xref->description;
        $description =~ s/'/`/g unless !defined $description;
        if (! defined $description) {
      $description = $pc_gene->description;
      if (defined $description) {
          $description =~ s/ \[Source.+//;
      }
      else {
          print STDERR "gene, " . $pc_gene->stable_id . ", without a description!\n";
      }
        }
        
        if ($description =~ /'/) {
      print SQL "INSERT INTO EG_Xref values (\\N,\"" . 
          $name . '",' .
          $tax_id. ',"' .
          $pc_gene->stable_id . '",' .
          $pc_gene->dbID . ',"' .
          $pcRNA->stable_id .'",' .
          $pcRNA->dbID .',"' .
          $pcRNA->biotype . '","' .
          $xref->database .'","' .
          $xref->display_id .'","' .
          $xref->primary_id .'","' .
          $description . '");'."\n" ;
        }
        else {
      print SQL "INSERT INTO EG_Xref values (\\N,'" . 
          $name . "'," .
          $tax_id. ",'" .
          $pc_gene->stable_id . "'," .
          $pc_gene->dbID . ",'" .
          $pcRNA->stable_id ."'," .
          $pcRNA->dbID .",'" .
          $pcRNA->biotype . "','" .
          $xref->database ."','" .
          $xref->display_id ."','" .
          $xref->primary_id ."','" .
          $description . "');\n" ;
        }
    }
      } # End foreach transcript
  } # End foreach gene
    }
    
    $dba->dbc->disconnect_if_idle;
    
}

1;
