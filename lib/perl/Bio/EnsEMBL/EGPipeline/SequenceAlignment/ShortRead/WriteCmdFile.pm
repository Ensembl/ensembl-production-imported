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

package Bio::EnsEMBL::EGPipeline::SequenceAlignment::ShortRead::WriteCmdFile;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Spec::Functions qw(catdir);
use JSON;

sub run {
  my ($self) = @_;
  
  my $species     = $self->param_required('species');
  my $aligner     = $self->param_required('aligner');
  my $results_dir = $self->param_required('results_dir');
  my $merge_level = $self->param_required('merge_level');
  my $merge_id    = $self->param_required('merge_id');
  my $bam_file    = $self->param_required('merged_bam_file');
  my $bw_file     = $self->param('bw_file');
  
  my %cmds;
  
  my $sql = "SELECT run_id, cmds, version FROM align_cmds WHERE merge_id = ? ORDER BY auto_id;";
  my $sth = $self->hive_dbh->prepare($sql);
  $sth->execute($merge_id);
  
  my %versions;
  my %run_ids;
  while (my $results = $sth->fetch) {
    my ($run_id, $cmds, $version) = @$results;
    
    $run_ids{$run_id}++ if $run_id;
    my @cmds = split(/\s*;\s*/, $cmds);
    foreach my $cmd (@cmds) {
      push @{$cmds{$merge_id}{'cmds'}}, $cmd;
      $versions{$version}++ if $version;
    }
  }
  push @{$cmds{$merge_id}{'run_ids'}}, sort keys %run_ids;
  
  # Get assembly name
  my $dba = $self->core_dba;
  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  
  $cmds{$merge_id}{'aligner'}         = $aligner;
  $cmds{$merge_id}{'aligner_version'} = join(", ", keys %versions);
  $cmds{$merge_id}{'merge_level'}     = $merge_level;
  $cmds{$merge_id}{'bam_file'}        = $bam_file;
  $cmds{$merge_id}{'bw_file'}         = $bw_file;
  
  my $json = to_json( \%cmds, { ascii => 1, pretty => 1 } );
  
  my $file_name = sprintf("/%s_%s.cmds.json", $merge_id, $assembly);
  my $cmds_file = catdir($results_dir, $file_name);
  open (my $fh, '>', $cmds_file) or die "Failed to open file '$cmds_file'";
	print $fh $json;
  close($fh);
  
  $self->param('cmds_file', $cmds_file);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id({'cmds_file' => $self->param('cmds_file')}, 1);
}

1;
