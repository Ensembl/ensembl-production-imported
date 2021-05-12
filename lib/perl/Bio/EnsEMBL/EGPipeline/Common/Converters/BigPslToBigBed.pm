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

package Bio::EnsEMBL::EGPipeline::Common::Converters::BigPslToBigBed;

use strict;
use warnings;
use Capture::Tiny ':all';
use JSON;

use File::Basename;
use File::Spec::Functions qw(catfile);

#use Bio::EnsEMBL::Utils::Exception qw(throw);

use base ('Bio::EnsEMBL::EGPipeline::FileDump::BaseDumper');

use constant _HINT_NAMES => [qw/type tab extraIndex/];

sub param_defaults {
  my ($self) = @_;
  
  return {
    'program'       => 'bedToBigBed',
    'file_varname'  => 'bigPsl_file', 
    'file_type'     => 'bb',
    'auto_sql_file' => catfile(dirname(__FILE__), 'bigPsl.as'),
    'bb_hints'      => [qw/ -type=bed12+13 -tab /],
  };
}

sub run {
  my ($self) = @_;
  
  my $program          = $self->param('program');
  my $chrom_sizes_file = $self->param_required('chrom_sizes_file');
  my $bigPsl_file         = $self->param_required('bigPsl_file');

  my $out_file         = $self->param('out_file');

  my $autosql_file     = $self->param('auto_sql_file');
  my $bb_hints         = $self->param('bb_hints');
  
  #First check bigPsl_file and chrom_sizes files
  if (! -s $bigPsl_file) {
    $self->warning("no or empty input bed file: $bigPsl_file" );
    return;
  }

  (-s $chrom_sizes_file)
    or $self->throw("no or empty input chrom sizes file: $chrom_sizes_file" );

  if (!defined $out_file) {
      ($out_file = $bigPsl_file) =~ s/\.bigPsl$/\.bb/;
  } 
  
  # sorting
  my $sorted_bigPsl_file = $bigPsl_file . '.sorted';
  my $sort_cmd = "LC_ALL=C sort -k1,1 -k2,2n $bigPsl_file > $sorted_bigPsl_file";
  $self->_execute($sort_cmd);

  # converting
  my @bb_cmd_args = ($program);
  #$self->push_hints(\@bb_cmd_args, $param_hints, $file_hints);
  push @bb_cmd_args, @$bb_hints; 
  push @bb_cmd_args, "-as='".$autosql_file."'" if (defined $autosql_file);
  push @bb_cmd_args, ($sorted_bigPsl_file, $chrom_sizes_file, $out_file);

  my $bb_cmd = join(' ', @bb_cmd_args);
  $self->warning("Trying to convert bigPsl to BigBed: $bb_cmd");
  $self->_execute($bb_cmd);

  $self->param('out_file', $out_file);
}

sub write_output {
  my ($self) = @_;
  
  my $dataflow_output = {
    'bigPsl_file' => $self->param('out_file'),
    'out_file' => $self->param('out_file'),
  };
  
  $self->dataflow_output_id($dataflow_output, 1);
}

sub _execute {
  my $self = shift;
  my ($cmd) = @_;
  
  
  my ($stdout, $stderr, $exit) = capture {
    system($cmd);
  };
  if ($exit) {
    $self->throw("Cannot execute $cmd:\n$stderr");
  }
}

1;
