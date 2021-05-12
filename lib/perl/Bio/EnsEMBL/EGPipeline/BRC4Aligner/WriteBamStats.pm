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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::WriteBamStats;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use Capture::Tiny ':all';
use Path::Tiny qw(path);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

sub write_output {
  my ($self) = @_;

  my $stats = $self->param_required('bam_stats');
  my $results_dir = $self->param_required('results_dir');

  my $stats_file = catdir($results_dir, 'mappingStats.txt');

  my @fields = qw(file coverage mapped number_reads_mapped average_read_length number_pairs_mapped);

  # Print header
  open my $statsh, ">", $stats_file;
  print $statsh join("\t", @fields) . "\n";
  for my $stat (sort { $a->{file} cmp $b->{file} } @$stats) {
    my @data_fields = map { $stat->{$_} } @fields;
    print $statsh join("\t", @data_fields) . "\n";
  }
  close $statsh;
}
1;
