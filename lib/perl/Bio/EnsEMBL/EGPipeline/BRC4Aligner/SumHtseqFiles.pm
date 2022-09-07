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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::SumHtseqFiles;
use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::BRC4Aligner::Base');

use Data::Dumper qw(Dumper);
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir catfile);


sub param_defaults {
  my ($self) = @_;

  return {
    # Features on which to count the read
    features => 'exon',
  };
}

sub run {
  my ($self) = @_;

  my $cases = $self->param('case');
  my @unstranded;
  my @uniquestranded;
  my @nonuniquestranded;
  my @sum_htseq_files;

  for my $case_hash (@$cases) {
    push @unstranded, $case_hash if $case_hash->{"strand"} eq 'unstranded';
    push @uniquestranded, $case_hash if ($case_hash->{"strand"} ne 'stranded') && ($case_hash->{"number"} eq 'unique');
    push @nonuniquestranded, $case_hash if ($case_hash->{"strand"} ne 'stranded') && ($case_hash->{"number"} eq 'total');
  }

  #Unstranded
  if (@unstranded != 0) {
    my @unstranded_htseq_files = map {$_->{htseq_file}} @unstranded;
    for my $unstranded_htseq_file (@unstranded_htseq_files) {
      push @sum_htseq_files, $unstranded_htseq_file
    };
  };

  if (@uniquestranded != 0) {
    my @uqstranded_htseq_files = map {$_->{htseq_file}} @uniquestranded;
    my $results_dir = dirname($uqstranded_htseq_files[0]);
    my $sum_uqstranded_htseq_file = catfile($results_dir, 'genes.htseq-union.stranded.sum.counts');
    my $uqstranded_htseq_files_str = join " ", @uqstranded_htseq_files;
    my $sum_uqstranded_htseq_file_output = `cat $uqstranded_htseq_files_str | awk '{a[\$1]+=\$2}END{for(i in a) print i,a[i]}' > $sum_uqstranded_htseq_file`;
    push @sum_htseq_files, $sum_uqstranded_htseq_file;
  };


  if (@nonuniquestranded != 0) {
    my @nuqstranded_htseq_files = map {$_->{htseq_file}} @nonuniquestranded;
    my $results_dir = dirname($nuqstranded_htseq_files[0]);
    my $sum_nuqstranded_htseq_file = catfile($results_dir, 'genes.htseq-union.stranded.nonunique.sum.counts');
    my $nuqstranded_htseq_files_str = join " ", @nuqstranded_htseq_files;
    my $sum_nuqstranded_htseq_file_output = `cat $nuqstranded_htseq_files_str | awk '{a[\$1]+=\$2}END{for(i in a) print i,a[i]}' > $sum_nuqstranded_htseq_file`;
    push @sum_htseq_files, $sum_nuqstranded_htseq_file;
  };

  for my $sum_htseq_file (@sum_htseq_files) {
    print($sum_htseq_file."\n");
  }

  $self->param("sum_htseq_files", \@sum_htseq_files);
}

sub write_output {

  my ($self) = @_;

  my @sum_htseq_files = @{ $self->param("sum_htseq_files") };

  my @cases;

  for my $sum_htseq_file (@sum_htseq_files) {
    push @cases, {
        sum_file => $sum_htseq_file,
    };
  }

  for my $case (@cases) {
    $self->dataflow_output_id($case, 2);
  }

}

1;