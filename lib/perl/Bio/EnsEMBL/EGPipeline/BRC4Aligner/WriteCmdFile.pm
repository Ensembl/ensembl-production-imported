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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::WriteCmdFile;

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
  my $sample_name    = $self->param_required('sample_name');
  
  
  my $sql = "SELECT cmds, version FROM align_cmds WHERE sample_name = ? ORDER BY auto_id;";
  my $sth = $self->hive_dbh->prepare($sql);
  $sth->execute($sample_name);
  
  my %cmds;
  my %versions;
  while (my $results = $sth->fetch) {
    my ($cmds, $version) = @$results;
    
    my @cmds = split(/\s*;\s*/, $cmds);
    foreach my $cmd (@cmds) {
      push @{ $cmds{'cmds'} }, $cmd;
      $versions{$version}++ if $version;
    }
  }
  
  # Get assembly name
  my $dba = $self->core_dba;
  my $assembly = $dba->get_MetaContainer()->single_value_by_key('assembly.default');
  
  $cmds{sample_name} = $sample_name;
  $cmds{'aligner'}         = $aligner;
  $cmds{'aligner_version'} = join(", ", keys %versions);
  
  my $json = to_json( \%cmds, { ascii => 1, pretty => 1 } );
  
  my $file_name = "commands.json";
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
