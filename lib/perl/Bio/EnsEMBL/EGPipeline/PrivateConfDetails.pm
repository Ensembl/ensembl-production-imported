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

PrivateConfDetails

=head1 DESCRIPTION

EG specific private configuration details.
A stub. There should be a loadable hash ref "$CONFIG" from 
  Bio::EnsEMBL::EGPipeline::PrivateConfDetails::Impl
for this to work properly.

=cut

package Bio::EnsEMBL::EGPipeline::PrivateConfDetails;

use strict;
use warnings;

our @EXPORT_OK = qw(private_conf);
use Exporter;
our @ISA = qw(Exporter);

use Try::Tiny;

try {
  require Bio::EnsEMBL::EGPipeline::PrivateConfDetails::Impl;
};

=head2 private_conf

Description: Method thats returns value from the private configuration for the specified "conf_key", or undef -- if not key met

=cut

sub private_conf {
  my ($conf_key) = @_;

  my $val = try {
    my $config = $Bio::EnsEMBL::EGPipeline::PrivateConfDetails::Impl::CONFIG;
    exists $config->{$conf_key} && $config->{$conf_key};
  } || undef;

  return $val;
}

1;
