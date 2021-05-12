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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::RepeatMaskerFactory;

use strict;
use warnings;

use Bio::SeqIO;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub param_defaults {
  return {
    'max_seq_length' => 1000000,
  };
}

sub write_output {
  my ($self) = @_;
  
  my $species = $self->param_required('species');
  my $always_use_repbase = $self->param_required('always_use_repbase');
  my $rm_library = $self->param_required('rm_library');
  my $rm_logic_name = $self->param_required('rm_logic_name');
  
  # Separate logic_names for separate repeat masker runs
  # By default, we use repbase, unless there is a custom rm library for this species
  my @logic_names;
  if ($always_use_repbase || (! exists $$rm_library{$species} && ! exists $$rm_library{'all'})) {
    push @logic_names, 'repeatmask_repbase';
  }

  # Define the logic_name for the custom rm library run
  if (exists $$rm_library{$species} || exists $$rm_library{'all'}) {
    if (exists $$rm_logic_name{$species} || exists $$rm_logic_name{'all'}) {
      my $name = $$rm_logic_name{$species} || $$rm_logic_name{'all'};
      push @logic_names, $name;
    } else {
      push @logic_names, 'repeatmask_customlib';
    }
  }
  
  my $queryfile = $self->param_required('queryfile');
  if (!-e $queryfile) {
    $self->throw("Query file '$queryfile' does not exist");
  } else {
    my $max_seq_length = $self->param('max_seq_length');
    
    # We want to split the query if it is too long
    my $split_up_query = 0;
    
    if (defined $max_seq_length) {
      my $total_length = 0;
      my $fasta = Bio::SeqIO->new(-format => 'Fasta', -file => $queryfile);

      while (my $seq = $fasta->next_seq) {
        $total_length += $seq->length;
      }

      if ($total_length > $max_seq_length) {
        $split_up_query = 1;
      }
    }

    if ($split_up_query) {
      my $dba = $self->get_DBAdaptor($self->param('db_type'));
      my $slice_adaptor = $dba->get_adaptor('Slice');

      my $fasta = Bio::SeqIO->new(-format => 'Fasta', -file => $queryfile);
      while (my $seq = $fasta->next_seq) {
        my $seq_length = $seq->length;
        my ($start, $end) = (1, $max_seq_length);
        $end = $seq_length if $end > $seq_length;

        while ($start <= $seq_length) {
          my $querylocation = $seq->id.":$start-$end";

          foreach my $logic_name (@logic_names) {
            $self->dataflow_output_id(
              {
                'logic_name'    => $logic_name,
                'queryfile'     => undef,
                'querylocation' => $querylocation,
              }, 1);
          }

          $start = $end + 1;
          $end = $start + $max_seq_length - 1;
          $end = $seq_length if $end > $seq_length;
        }
      }
      $dba->dbc->disconnect_if_idle();

    } else {
      foreach my $logic_name (@logic_names) {
        $self->dataflow_output_id(
          {
            'logic_name' => $logic_name,
            'queryfile'  => $queryfile,
          }, 1);
      }
    }
  }
}

1;
