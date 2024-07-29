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

package Bio::EnsEMBL::EGPipeline::DNAFeatures::Utils;

use strict;
use warnings;
use File::Spec::Functions qw(catfile catdir);
use File::Basename qw(basename dirname);

our @EXPORT_OK = qw(get_repbase_species_lib);
use Exporter;
our @ISA = qw(Exporter);

sub get_repbase_species_lib {
  my ($dba, $lib_path) = @_;
  
  # Get Repeat masker script path to run
  my $script = _get_repeat_script_path();
  print STDERR "RepeatMasker script is at $script\n";
  
  # Get list of taxonomy rank for the given species
  my $classification = _get_classification($dba);
  print STDERR "Classification = " . join(", ", @$classification) . "\n";

  my $temp_lib_path = $lib_path . ".tmp";
  
  my $best_tax_name = '';
  for my $tax_name (@$classification) {
    my $command = "perl $script -species '$tax_name' > $temp_lib_path";
    system($command);
    
    if (-s $temp_lib_path) {
      $best_tax_name = $tax_name;
      print STDERR "Repbase taxonomy used: $tax_name\n";
      last;
    }
  }
  die "Could not generate a Repbase lib for $classification->[0]" if not -s $temp_lib_path;
  
  rename $temp_lib_path, $lib_path;
  
  # Log the tax name used
  my $lib_species_path = $lib_path . ".species";
  open my $lib_sp, ">", $lib_species_path;
  print $lib_sp "RepBase taxonomy: $best_tax_name\n";
  close $lib_sp;
  
  return $best_tax_name;
}

sub _get_classification {
  my ($dba) = @_;
  
  my $ma = $dba->get_adaptor('MetaContainer');

  my ($scientific_name) = @{ $ma->list_value_by_key('species.scientific_name') };
  my @classification = @{ $ma->list_value_by_key('species.classification') };
  my ($genus) = split " ", $scientific_name;
  
  unshift @classification, ($scientific_name, $genus);
  
  print "CLASSIFICATION: \n" . join("\n", @classification) . "\n";
  
  return \@classification;
}

sub _get_repeat_script_path {
#  my () = @_;
  
  my $rm_path = `which RepeatMasker`;
  chomp $rm_path;
  die "Could not find RepeatMasker in your path" unless $rm_path;
  
  # Resolve symlink
  if (-l $rm_path) {
    my $target_path = readlink $rm_path;
    if ($target_path =~ /^\//) {
      $rm_path = $target_path;
    } else {
      $rm_path = catfile(dirname($rm_path), $target_path);
    }
  }
  
  # Get to the util dir
  my $rm_dir = dirname($rm_path);

  my $util_dir = catdir($rm_dir, '../libexec/util');
  $util_dir = catdir($rm_dir, '../lib/util') if not -d $util_dir;
  die "Could not find RepeatMasker util dir at $util_dir" if not -d $util_dir;
  
  my $script_path = catfile($util_dir, 'queryRepeatDatabase.pl');
  die "Could not find RepeatMasker queryRepeatDatabase.pl script" if not -e $script_path;
  
  return $script_path;
}

1;
