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

package Bio::EnsEMBL::EGPipeline::Common::RunnableDB::DumpProteome;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Bio::EnsEMBL::Utils::IO::FASTASerializer;

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'overwrite'         => 0,
    'header_style'      => 'default',
    'chunk_factor'      => 1000,
    'line_width'        => 80,
    'allow_stop_codons' => 0,
    'is_canonical'      => undef,
    'production_lookup' => 0,
    'file_varname'      => 'proteome_file',
  };
  
}

sub fetch_input {
  my ($self) = @_;
  
  my $proteome_file = $self->param('proteome_file');
  my $proteome_dir  = $self->param('proteome_dir');
  my $species       = $self->param('species');
  
  if (!defined $proteome_file) {
    if (!defined $proteome_dir) {
      $self->throw("A path or filename is required");
    } else {
      if (!-e $proteome_dir) {
        $self->warning("Output directory '$proteome_dir' does not exist. I shall create it.");
        make_path($proteome_dir) or $self->throw("Failed to create output directory '$proteome_dir'");
      }
      $proteome_file = catdir($proteome_dir, "$species.fa");
      $self->param('proteome_file', $proteome_file);
    }
  }
  
  if (-e $proteome_file) {
    # TODO: Could count the number of header lines here and check against the number of
    # protein_coding transcripts to help catch cases where the file exists but is incomplete or empty.
    if ($self->param('overwrite')) {
      $self->warning("Proteome file '$proteome_file' already exists, and will be overwritten.");
    } else {
      $self->warning("Proteome file '$proteome_file' already exists, and won't be overwritten.");
      $self->param('skip_dump', 1);
    }
  }
}

sub run {
  my ($self) = @_;
  
  return if $self->param('skip_dump');
  
  my $proteome_file     = $self->param('proteome_file');
  my $header_style      = $self->param('header_style');
  my $chunk_factor      = $self->param('chunk_factor');
  my $line_width        = $self->param('line_width');
  my $allow_stop_codons = $self->param('allow_stop_codons');
  my $is_canonical      = $self->param('is_canonical');
  my $production_lookup = $self->param('production_lookup');
  
  # Use the ensembl_production database to retrieve the biotypes
  # associated with the coding group, if required.
  my $biotypes;
  
  if ($production_lookup) {
    my $biotype_groups = ['coding'];
    my $pdba = $self->get_DBAdaptor('production');
    my $biotype_manager = $pdba->get_biotype_manager();
    map { push @{$biotypes}, @{ $biotype_manager->group_members($_)} } @{$biotype_groups};
  } else {
    push @{$biotypes}, 'protein_coding';
  }
  
  open(my $fh, '>', $proteome_file) or $self->throw("Cannot open file $proteome_file: $!");
  my $serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new(
    $fh,
    undef,
    $chunk_factor,
    $line_width,
  );
  
  my $dba = $self->core_dba();
  my $tra = $dba->get_adaptor('Transcript');
  my $transcripts = $tra->fetch_all_by_biotype($biotypes);
  
  foreach my $transcript (sort { $a->seq_region_name cmp $b->seq_region_name or $a->seq_region_start <=> $b->seq_region_start } @{$transcripts}) {
    if (defined $is_canonical) {
      next if $is_canonical != $transcript->is_canonical;
    }
    
    my $seq_obj = $transcript->translate();
    
    if ($header_style ne 'default') {
      $seq_obj->display_id($self->header($header_style, $transcript));
    }
    
    if ($seq_obj->seq() =~ /\*/ && !$allow_stop_codons) {
      $self->warning("Translation for transcript ".$transcript->stable_id." contains stop codons. Skipping.");
    } else {
      if ($seq_obj->seq() =~ /\*/) {
        $self->warning("Translation for transcript ".$transcript->stable_id." contains stop codons.");
      }
      $serializer->print_Seq($seq_obj);
    }
	}
  $dba->dbc->disconnect_if_idle();
  
  close($fh);
}

sub write_output {
  my ($self) = @_;
  my $file_varname = $self->param_required('file_varname');
  
  $self->dataflow_output_id({$file_varname => $self->param('proteome_file')}, 1);
}

sub header {
  my ($self, $header_style, $transcript) = @_;
  
  my $translation = $transcript->translation;
  my $header = $translation->stable_id;
  
  if ($header_style eq 'dbID') {
    $header = $translation->dbID;
    
  } elsif ($header_style eq 'extended') {
    my $gene = $transcript->get_Gene;
    my $id = $translation->stable_id;
    my $desc = $gene->description ? $gene->description : ' ';
    $desc =~ s/\s*\[Source.+$//;
    
    my $xref = $gene->display_xref;
    if (defined $xref) {
      my $name = $xref->display_id;
      $desc = "$name: $desc";
    }
    
    my $location = join(':',
      $transcript->seq_region_name,
      $translation->genomic_start,
      $translation->genomic_end,
      $transcript->strand,
    );
    
    $header = join('|',
      "$id $desc",
      $transcript->biotype,
      $location,
      'gene:'.$gene->stable_id
    );
    
  }
  
  return $header;
}

1;
