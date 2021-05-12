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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::BlastEG

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Runnable::BlastEG->
  new(
      -query    => $slice,
      -program  => 'blastn',
      -database => 'uniprot_sprot_invertebrates',
      -options  => 'hitdist=40',
      -parser   => $parser,
      -filter   => undef,
     );
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

  This module is a wrapper for running ncbi-blast. 
  It requires a parser object, e.g. FilterBPlite.

=cut

package Bio::EnsEMBL::Analysis::Runnable::BlastEG;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($database, $parser) = rearrange(['DATABASE', 'PARSER'], @args);
  
  throw("DATABASE option is mandatory") unless defined $database;
  throw("PARSER option is mandatory") unless defined $parser;
  
  $self->database($database);
  $self->parser($parser);
  
  return $self;
}

sub database {
  my $self = shift;
  $self->{'database'} = shift if(@_);
  return $self->{'database'};
}

sub parser {
  my $self = shift;
  $self->{'parser'} = shift if(@_);
  return $self->{'parser'};
}

sub run_analysis {
  my ($self) = @_;
  my $command     = $self->program;
  my $database    = $self->database;
  my $queryfile   = $self->queryfile;
  my $resultsfile = $self->resultsfile;

  $command .= " -db $database -query $queryfile -out $resultsfile ";
  $command .= $self->options;
  
  system($command) == 0 or throw("FAILED to run ".$command);
}

sub parse_results {
  my ($self, $resultsfile) = @_;
  $resultsfile = $self->resultsfile unless defined $resultsfile;
  
  # Remove NCBI footer as it confuses the Ensembl parser.
  my $sed_command = "sed -i '/^Lambda/, \$d' $resultsfile";
  system($sed_command) == 0 or throw("FAILED to run ".$sed_command);
  
  $self->output($self->parser->parse_files([$resultsfile]));
}

1;
