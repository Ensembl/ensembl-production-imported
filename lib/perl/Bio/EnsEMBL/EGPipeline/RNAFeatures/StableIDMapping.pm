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

=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::RNAFeatures::StableIDMapping

=cut

package Bio::EnsEMBL::EGPipeline::RNAFeatures::StableIDMapping;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);
use POSIX qw(strftime);

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  my ($self) = @_;
  return {
    'report_style' => 'eg',
  };
}

sub fetch_input {
  my ($self) = @_;
  my $species      = $self->param_required('species');
  my $report_style = $self->param_required('report_style');
  my $report_dir   = $self->param_required('report_dir');
  
  my $old_db = {
    -host    => $self->param_required('old_host'),
    -port    => $self->param_required('old_port'),
    -user    => $self->param_required('old_user'),
    -pass    => $self->param('old_pass') || '',
    -dbname  => $self->param_required('old_dbname'),
    -species => $species,
  };
  my $old_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$old_db);
  $self->param('old_dba', $old_dba);
  
  my $new_dba = $self->core_dba();
  $self->param('new_dba', $new_dba);
  
  my $old_mca = $old_dba->get_adaptor('MetaContainer');
  my $old_genebuild = $old_mca->single_value_by_key('genebuild.version');
  $self->param('old_genebuild', $old_genebuild);
  
  my $new_mca = $new_dba->get_adaptor('MetaContainer');
  my $new_genebuild = $new_mca->single_value_by_key('genebuild.version');
  $self->param('new_genebuild', $new_genebuild);
  
  my $report_file;
  if ($report_style eq 'vb') {
    my $vb_species = $new_mca->single_value_by_key('species.vectorbase_name');
    my $strain     = $new_mca->single_value_by_key('species.strain');
    $report_file   = ucfirst("$vb_species-$strain");
  } else {
    $report_file = $species;
  }
  $report_file .= "_MAPPINGS_$old_genebuild-$new_genebuild.txt";
  
  $self->param('mapping_file', catdir($report_dir, $report_file));
}

sub run {
  my ($self) = @_;
  
  $self->generate_mapping_session();
  $self->find_mappings();
  $self->save_mapping_report();
}

sub generate_mapping_session {
  my ($self) = @_;
  my $old_dba = $self->param_required('old_dba');
  my $new_dba = $self->param_required('new_dba');
  
  my $max_id_sql = 'SELECT MAX(mapping_session_id) FROM mapping_session';
  my $sth = $new_dba->dbc->prepare($max_id_sql) or $self->throw($new_dba->dbc->db_handle->errstr);
  $sth->execute() or $self->throw($sth->errstr);
  my ($mapping_session_id) = $sth->fetchrow_array();
  $mapping_session_id++;
  
  my $old_dbname    = $old_dba->dbc->dbname();
  my $new_dbname    = $new_dba->dbc->dbname();
  my ($old_release) = $old_dbname =~ /core_(\d+)/;
  my ($new_release) = $new_dbname =~ /core_(\d+)/;
  my $old_genebuild = $self->param_required('old_genebuild');
  my $new_genebuild = $self->param_required('new_genebuild');
  
  my $insert_sql = 
    'INSERT INTO mapping_session ('.
      'mapping_session_id, '.
      'old_db_name, '.
      'new_db_name, '.
      'old_release, '.
      'new_release, '.
      'old_assembly, '.
      'new_assembly, '.
      'created) '. 
    'VALUES (?, ?, ?, ?, ?, ?, ?, ?)';
  $sth = $new_dba->dbc->prepare($insert_sql) or $self->throw($new_dba->dbc->db_handle->errstr);
  $sth->execute(
    $mapping_session_id,
    $old_dbname,
    $new_dbname,
    $old_release,
    $new_release,
    $old_genebuild,
    $new_genebuild,
    strftime("%Y-%m-%d %T", localtime(time))
  ) or $self->throw($sth->errstr);
  
  $self->param('mapping_session_id', $mapping_session_id);
}

sub find_mappings {
  my ($self, ) = @_;
  my $old_dba = $self->param_required('old_dba');
  my $new_dba = $self->param_required('new_dba');
  
  my $old_ta  = $old_dba->get_adaptor('Transcript');
  my $new_ta  = $new_dba->get_adaptor('Transcript');
  
  # Mapping is somewhat restricted in this case, so we don't have to
  # consider all possible scenarios. There will be a one-to-one mapping
  # between gene and transcript. If a stable ID is the same in new and
  # old dbs, then it hasn't changed (the ID generating code will always
  # create a new ID if there is not an exact match to an existing ID).
  # So we want to get two sets of stable IDs, those only in the old db
  # and those only in the new db. Then we'll look at the overlap
  # between the sets and do mapping accordingly.
  my %old_transcripts;
  foreach my $transcript (@{$old_ta->fetch_all}) {
    $old_transcripts{$transcript->stable_id} = $transcript;
  }
  
  my %new_transcripts;
  foreach my $transcript (@{$new_ta->fetch_all}) {
    if (exists $old_transcripts{$transcript->stable_id}) {
      delete $old_transcripts{$transcript->stable_id};
    } else {
      $new_transcripts{$transcript->stable_id} = $transcript;
    }
  }
  
  my %seen;
  foreach my $old_id (sort keys %old_transcripts) {
    my $mapping = 0;  
    my $overlaps = $new_ta->fetch_all_nearest_by_Feature(
      -FEATURE => $old_transcripts{$old_id},
      -RANGE => 0,
    );
    foreach my $overlap (@$overlaps) {
      my $new_transcript = $$overlap[0];
      if (exists $new_transcripts{$new_transcript->stable_id}) {
        $self->add_mapping($old_transcripts{$old_id}, $new_transcript);
        $seen{$new_transcript->stable_id}++;
        $mapping = 1;
      }
    }
    if (!$mapping) {
      $self->add_mapping($old_transcripts{$old_id}, undef);
    }
  }
  
  foreach my $new_id (sort keys %new_transcripts) {
    if (! exists $seen{$new_id}) {
      $self->add_mapping(undef, $new_transcripts{$new_id});
    }
  }
}

sub add_mapping {
  my ($self, $old_transcript, $new_transcript) = @_;
  
  my ($old_t_stable_id, $old_t_version, $old_g_stable_id, $old_g_version);
  my ($new_t_stable_id, $new_t_version, $new_g_stable_id, $new_g_version);
  
  if (defined $old_transcript) {
    my $old_gene     = $old_transcript->get_Gene();
    $old_t_stable_id = $old_transcript->stable_id;
    $old_t_version   = $old_transcript->version;
    $old_g_stable_id = $old_gene->stable_id;
    $old_g_version   = $old_gene->version;
  }
  
  if (defined $new_transcript) {
    my $new_gene     = $new_transcript->get_Gene();
    $new_t_stable_id = $new_transcript->stable_id;
    $new_t_version   = $new_transcript->version;
    $new_g_stable_id = $new_gene->stable_id;
    $new_g_version   = $new_gene->version;
  }
  
  if (defined $old_transcript) {    
    $self->add_gene_archive(
      $old_g_stable_id, $old_g_version, $old_t_stable_id, $old_t_version);
    
    if (defined $new_transcript) {
      $self->add_stable_id_event('Transcript', $old_t_stable_id, $old_t_version, $new_t_stable_id, $new_t_version);
      $self->add_stable_id_event('Gene', $old_g_stable_id, $old_g_version, $new_g_stable_id, $new_g_version);
    } else {
      $self->add_stable_id_event('Transcript', $old_t_stable_id, $old_t_version, undef, undef);
      $self->add_stable_id_event('Gene', $old_g_stable_id, $old_g_version, undef, undef);
    }
  } elsif (defined $new_transcript) {
    $self->add_stable_id_event('Transcript', undef, undef, $new_t_stable_id, $new_t_version);
    $self->add_stable_id_event('Gene', undef, undef, $new_g_stable_id, $new_g_version);
  }
  
  $self->mapping_report_line($old_transcript, $new_transcript);
}

sub add_stable_id_event {
  my ($self, $type, $old_stable_id, $old_version, $new_stable_id, $new_version) = @_;
  
  my $new_dba            = $self->param_required('new_dba');
  my $mapping_session_id = $self->param_required('mapping_session_id');
  
  my $sql =
    'INSERT IGNORE INTO stable_id_event ('.
      'mapping_session_id, '.
      'type, '.
      'old_stable_id, '.
      'old_version, '.
      'new_stable_id, '.
      'new_version) '.
    'VALUES (?, ?, ?, ?, ?, ?);';
  
  my $sth = $new_dba->dbc->prepare($sql) or $self->throw($new_dba->dbc->db_handle->errstr);
  $sth->execute(
    $mapping_session_id,
    $type,
    $old_stable_id,
    $old_version,
    $new_stable_id,
    $new_version
  ) or $self->throw($sth->errstr);
}

sub add_gene_archive {
  my ($self, $old_g_stable_id, $old_g_version, $old_t_stable_id, $old_t_version) = @_;
  
  my $new_dba            = $self->param_required('new_dba');
  my $mapping_session_id = $self->param_required('mapping_session_id');
  
  my $sql =
    'INSERT INTO gene_archive ('.
      'mapping_session_id, '.
      'gene_stable_id, '.
      'gene_version, '.
      'transcript_stable_id, '.
      'transcript_version) '.
    'VALUES (?, ?, ?, ?, ?);';
  
  my $sth = $new_dba->dbc->prepare($sql) or $self->throw($new_dba->dbc->db_handle->errstr);
  $sth->execute(
    $mapping_session_id,
    $old_g_stable_id,
    $old_g_version || 1,
    $old_t_stable_id,
    $old_t_version || 1
  ) or $self->throw($sth->errstr);
}

sub mapping_report_header {
  my ($self) = @_;
  
  my @header = qw(
    old_stable_id
    old_biotype
    old_seq_region
    old_start
    old_end
    old_strand
    new_stable_id
    new_biotype
    new_seq_region
    new_start
    new_end
    new_strand
  );
  
  return join("\t", @header)."\n";
}

sub mapping_report_line {
  my ($self, $old_transcript, $new_transcript) = @_;
  
  my $report = $self->param('mapping_report');
  if (! defined $report) {
    $report = $self->mapping_report_header();
  }
  
  my @old;
  my @new;
  
  if (defined $old_transcript) {
    @old = (
      $old_transcript->stable_id,
      $old_transcript->biotype,
      $old_transcript->seq_region_name,
      $old_transcript->seq_region_start,
      $old_transcript->seq_region_end,
      $old_transcript->seq_region_strand,
    );
  } else {
    @old = ('') x 6;
  }
  
  if (defined $new_transcript) {
    @new = (
      $new_transcript->stable_id,
      $new_transcript->biotype,
      $new_transcript->seq_region_name,
      $new_transcript->seq_region_start,
      $new_transcript->seq_region_end,
      $new_transcript->seq_region_strand,
    );
  } else {
    @new = ('') x 6;
  }
  
  $report .= join("\t", @old, @new)."\n";
  
  $self->param('mapping_report', $report);
}

sub save_mapping_report {
  my ($self) = @_;
  my $mapping_file   = $self->param_required('mapping_file');
  my $mapping_report = $self->param_required('mapping_report');
  
  my $file = path($mapping_file);
  $file->spew($mapping_report);
}

1;
