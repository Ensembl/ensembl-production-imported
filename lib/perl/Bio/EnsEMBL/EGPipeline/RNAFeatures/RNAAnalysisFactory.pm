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

package Bio::EnsEMBL::EGPipeline::RNAFeatures::RNAAnalysisFactory;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $species             = $self->param_required('species');
  my $analyses            = $self->param_required('analyses');
  my $run_cmscan          = $self->param_required('run_cmscan');
  my $run_trnascan        = $self->param_required('run_trnascan');
  my $load_mirbase        = $self->param_required('load_mirbase');
  my $rfam_logic_name     = $self->param_required('rfam_logic_name');
  my $cmscan_cm_file      = $self->param_required('cmscan_cm_file');
  my $cmscan_logic_name   = $self->param_required('cmscan_logic_name');
  my $trnascan_logic_name = $self->param_required('trnascan_logic_name');
  my $mirbase_logic_name  = $self->param_required('mirbase_logic_name');
  my $mirbase_files       = $self->param_required('mirbase_files');
  my $pipeline_dir        = $self->param_required('pipeline_dir');
  my $db_backup_file      = $self->param('db_backup_file');
  
  my $filtered_analyses = [];
  foreach my $analysis (@{$analyses}) {
    my $logic_name = $$analysis{'logic_name'};
    
    if ($run_cmscan) {
      if ($logic_name eq $rfam_logic_name) {
        # If we don't have a given CM file, we'll be using the Rfam default.
        if (! exists $$cmscan_cm_file{$species} && ! exists $$cmscan_cm_file{'all'}) {
          push @$filtered_analyses, $analysis;
        }
        
      } elsif ($logic_name eq 'cmscan_custom') {
        if (exists $$cmscan_cm_file{$species} || exists $$cmscan_cm_file{'all'}) {
          my $cm_file = $$cmscan_cm_file{$species} || $$cmscan_cm_file{'all'};
          $$analysis{'db_file'} = $cm_file;
          
          if (exists $$cmscan_logic_name{$species} || exists $$cmscan_logic_name{'all'}) {
            my $name = $$cmscan_logic_name{$species} || $$cmscan_logic_name{'all'};
            $$analysis{'logic_name'} = $name;
          }
          
          push @$filtered_analyses, $analysis;
        }
      }
      
    }
    
    if ($run_trnascan) {
      if ($logic_name eq $trnascan_logic_name) {
        push @$filtered_analyses, $analysis;
      }
    }
    
    if ($load_mirbase) {
      if ($logic_name eq $mirbase_logic_name) {
        if (exists $$mirbase_files{$species}) {
          $$analysis{'db_file'} = $$mirbase_files{$species};
          push @$filtered_analyses, $analysis;
        }
      }
    }
  }
  
  foreach my $analysis (@{$filtered_analyses}) {
    $$analysis{'species'} = $species;
    $$analysis{'db_backup_file'} = $db_backup_file if defined $db_backup_file;
  }
  
  $self->param('filtered_analyses', $filtered_analyses);
}

sub write_output {
  my ($self) = @_;
  
  $self->dataflow_output_id($self->param('filtered_analyses'), 2);
}

1;
