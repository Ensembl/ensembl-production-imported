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

package Bio::EnsEMBL::EGPipeline::BRC4Aligner::CleanupBam;

use strict;
use warnings;
use Cwd;

use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my ($self) = @_;
  
  my $skip_cleanup = $self->param("skip_cleanup");
  return if $skip_cleanup;
  
  my $dir = $self->param("results_dir");
  my $pwd = getcwd();
  chdir $dir or die("Could not change to $dir from $pwd");
  
  my @files = <*.bam*>;
  
  my $nfiles = scalar @files;
  die("There are no bam files to delete") if $nfiles == 0;
  print("$nfiles bam files in $dir\n");
  
  for my $file (@files) {
    next if $file eq "results_name.bam";
    print "TO DELETE: $file\n";
    unlink($file) or die "Could not delete $file: $!";
  }
  
  chdir $pwd;
}

1;
