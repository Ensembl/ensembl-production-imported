#!/usr/bin/env/perl
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

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::PredictionExon;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Data::Dumper;

# Set up mapping session id (add it manually to mapping_session)
my $map_session_id = 1;

my ($from_host, $from_port, $from_user, $from_pass, $from_dbname, $from_species);
my ($host, $port, $user, $pass, $dbname, $species);

&GetOptions(
  'from_host=s'      => \$from_host,
  'from_port=s'      => \$from_port,
  'from_user=s'      => \$from_user,
  'from_pass=s'      => \$from_pass,
  'from_dbname=s'    => \$from_dbname,
  'from_species=s'   => \$from_species,

  'host=s'      => \$host,
  'port=s'      => \$port,
  'user=s'      => \$user,
  'pass=s'      => \$pass,
  'dbname=s'    => \$dbname,
  'species=s'   => \$species,
);

# Old release core database
my $old_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host    => $from_host,
    -port    => $from_port,
    -user    => $from_user,
    -pass    => $from_pass,
    -dbname  => $from_dbname,
    -species => $from_species,
);

#  New release core database with same assembly (or assembly with same seq_regions)
my $new_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host    => $host,
    -port    => $port,
    -user    => $user,
    -pass    => $pass,
    -dbname  => $dbname,
    -species => $species,
);

# Check if the toplevel sequence that have the same name in both databases have the same length

sub load_slices {
  my ($dba) = @_;
  my %slices;

  my $slice_adaptor = $dba->get_adaptor("Slice");
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
    my $seq_region_name = $slice->seq_region_name;
    $slices{$seq_region_name} = $slice;
  }
  return %slices;
}

my %new_slices = load_slices($new_dba);
my %old_slices = load_slices($old_dba);

foreach my $new_slice_name ( keys %new_slices ) {
    if ( exists $old_slices{$new_slice_name}) {
        if ($old_slices{$new_slice_name}->length == $new_slices{$new_slice_name}->length) {
        } else {
            print $new_slice_name . " NO MATCH\n";
            exit;
        }
    }
}

# Get genebuild versions

my $new_genebuild = 2;
my $old_genebuild = 1;

print STDERR "Old database genebuild: " . $old_genebuild . "\n";
print STDERR "New database genebuild: " . $new_genebuild . "\n";


# Get all the genes for the new database

my $type = "Gene";

my $new_gene_adaptor = $new_dba->get_adaptor("Gene");
my $old_gene_adaptor = $old_dba->get_adaptor("Gene");
$type = lc($type);
my @new_genes = @{$new_gene_adaptor->fetch_all};
my @old_genes = @{$old_gene_adaptor->fetch_all};
my @new_genes_stable_id_left = map {$_->stable_id} @new_genes;
my @old_genes_stable_id_left = map {$_->stable_id} @old_genes;

print STDERR "Number of genes in old database: " . scalar(@old_genes) . "\n";
print STDERR "Number of genes in new database: " . scalar(@new_genes) . "\n";

my $new_transcript_adaptor = $new_dba->get_adaptor("Transcript");
my $old_transcript_adaptor = $old_dba->get_adaptor("Transcript");
my @new_transcripts = @{$new_transcript_adaptor->fetch_all};
my @old_transcripts = @{$old_transcript_adaptor->fetch_all};
my @new_transcripts_stable_id_left = map {$_->stable_id} @new_transcripts;
my @old_transcripts_stable_id_left = map {$_->stable_id} @old_transcripts;

print STDERR "Number of transcripts in old database: " . scalar(@old_transcripts) . "\n";
print STDERR "Number of transcripts in new database: " . scalar(@new_transcripts) . "\n";

my $new_translation_adaptor = $new_dba->get_adaptor("Translation");
my $old_translation_adaptor = $old_dba->get_adaptor("Translation");
#print "Got adaptor\n";
my @new_translations_stable_id_left = @{$new_translation_adaptor->list_stable_ids()};
my @old_translations_stable_id_left = @{$old_translation_adaptor->list_stable_ids()};
#print "Fetched all\n";

print STDERR "Number of translations in old database: " . scalar(@old_translations_stable_id_left) . "\n";
print STDERR "Number of translations in new database: " . scalar(@new_translations_stable_id_left) . "\n";

print "INSERT INTO stable_id_event (old_stable_id,old_version,new_stable_id,new_version,mapping_session_id,type,score) VALUES ";

for my $gene (@new_genes) {
    
    my $projected_gene = $gene->project('chromosome', 'GCA_000150675.1');

    if ($projected_gene->[0]) {

        my $clone = $projected_gene->[0]->to_Slice();

        my $new_gene_stable_id = $gene->stable_id;
        my $seq_region_name =  $clone->seq_region_name;
        my $seq_region_strand = $clone->strand;
        my $seq_region_start = $clone->start;
        my $seq_region_end = $clone->end;

        my $old_slice_adaptor = $old_dba->get_adaptor("Slice");
        my $old_slice = $old_slice_adaptor->fetch_by_region( 'toplevel', $seq_region_name, $seq_region_start, $seq_region_end );

        my @old_genes = @{$old_slice->get_all_Genes};

        if (scalar (@old_genes) == 1) {
            my $old_gene_stable_id = $old_genes[0]->stable_id;
            my $old_gene = $old_genes[0];
            @new_genes_stable_id_left = grep { $_ ne $new_gene_stable_id } @new_genes_stable_id_left;
            @old_genes_stable_id_left = grep { $_ ne $old_gene_stable_id } @old_genes_stable_id_left;
            print "('$old_gene_stable_id',$old_genebuild,'$new_gene_stable_id',$new_genebuild,$map_session_id,'gene',1),";
            my @new_transcripts = @{$gene->get_all_Transcripts()};
            my @old_transcripts = @{$old_gene->get_all_Transcripts()};
            for my $transcript (@new_transcripts) {
                for my $old_transcript (@old_transcripts) {
                    #print STDERR $transcript->translateable_seq() . "\n";
                    #print STDERR $old_transcript->translateable_seq() . "\n";
                    if ($old_transcript->translateable_seq() && $transcript->translateable_seq()) {
                        my $new_transcript_stable_id = $transcript->stable_id;
                        my $old_transcript_stable_id = $old_transcript->stable_id;
                        @new_transcripts_stable_id_left = grep { $_ ne $new_transcript_stable_id } @new_transcripts_stable_id_left;
                        @old_transcripts_stable_id_left = grep { $_ ne $old_transcript_stable_id } @old_transcripts_stable_id_left;
                        print "('$old_transcript_stable_id',$old_genebuild,'$new_transcript_stable_id',$new_genebuild,$map_session_id,'transcript',1),";
                        if ( $transcript->translation && $old_transcript->translation) {
                            my $new_translation_stable_id = $transcript->translation->stable_id;
                            my $old_translation_stable_id = $old_transcript->translation->stable_id;
                            @new_translations_stable_id_left = grep { $_ ne $new_translation_stable_id } @new_translations_stable_id_left;
                            @old_translations_stable_id_left = grep { $_ ne $old_translation_stable_id } @old_translations_stable_id_left;
                            print "('$old_translation_stable_id',$old_genebuild,'$new_translation_stable_id',$new_genebuild,$map_session_id,'translation',1),";
                        }
                    }
                    next;
                }
            }
        }
    }
}


my $string = "";

for my $stable_id (@new_genes_stable_id_left) {
    $string = $string . "(NULL,$old_genebuild,'$stable_id',$new_genebuild,$map_session_id,'gene',1),";
}
for my $stable_id (@old_genes_stable_id_left) {
    $string = $string . "('$stable_id',$old_genebuild,NULL,$new_genebuild,$map_session_id,'gene',1),";
}
for my $stable_id (@new_transcripts_stable_id_left) {
    $string = $string . "(NULL,$old_genebuild,'$stable_id',$new_genebuild,$map_session_id,'transcript',1),";
}
for my $stable_id (@old_transcripts_stable_id_left) {
    $string = $string . "('$stable_id',$old_genebuild,NULL,$new_genebuild,$map_session_id,'transcript',1),";
}
for my $stable_id (@new_translations_stable_id_left) {
    $string = $string . "(NULL,$old_genebuild,'$stable_id',$new_genebuild,$map_session_id,'translation',1),";
}
for my $stable_id (@old_translations_stable_id_left) {
    $string = $string . "('$stable_id',$old_genebuild,NULL,$new_genebuild,$map_session_id,'translation',1),";
}

chop $string;

print $string . ";\n";

print STDERR "Genes removed in old db " . scalar(@old_genes_stable_id_left) . "\n";
print STDERR "Unmapped genes in new db " . scalar(@new_genes_stable_id_left) . "\n";
print STDERR "Transcripts removed in old db " . scalar(@old_transcripts_stable_id_left) . "\n";
print STDERR "Unmapped transcripts in new db " . scalar(@new_transcripts_stable_id_left) . "\n";
print STDERR "Transcripts removed in old db " . scalar(@old_translations_stable_id_left) . "\n";
print STDERR "Unmapped transcripts in new db " . scalar(@new_translations_stable_id_left) . "\n";


# For each gene get it's location, use that same location to find the corresponding gene in the old database
