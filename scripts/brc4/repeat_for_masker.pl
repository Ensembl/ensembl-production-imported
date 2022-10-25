#!/usr/env perl
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use v5.14.00;
use strict;
use warnings;
use Carp;
use autodie qw(:all);
use Readonly;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl qw( :easy ); 
Log::Log4perl->easy_init($WARN); 
my $logger = get_logger(); 

use Bio::EnsEMBL::Registry;
use JSON;
use File::Spec::Functions qw(rel2abs catfile);

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $ref_libs = get_libs($opt{libs_dir});
my $tab = get_tab_file($opt{tab});
my $ref_map = get_species_map($opt{ref_registry});
my $new_map = get_species_map($opt{new_registry});

my $libs = get_libraries($tab, $ref_libs, $ref_map, $new_map);

# Get the list of species to mask before removing the empty custom libs
my @species = sort keys %$libs;
$libs = remove_empty_libs($libs);

my $output = format_json_for_masker(\@species, $libs);

write_json($output, $opt{libs_dir}, $opt{out_json});

###############################################################################

sub format_json_for_masker {
  my ($species, $libs) = @_;
  
  my %data = (
    species => $species,
    repeatmasker_library => $libs,
  );
  
  return \%data;
}

sub remove_empty_libs {
  my ($libs) = @_;
  
  for my $sp (keys %$libs) {
    my $file = $libs->{$sp};
    if (not -s $file) {
      warn("Empty repeat library for $sp (no custom lib will be used): $file\n");
      delete $libs->{$sp};
    }
  }
  return $libs;
}

sub get_libraries {
  my ($tab, $ref, $map, $to_mask) = @_;

  my %libraries;
  for my $entry (@$tab) {
    # Get the species
    my $organism = $entry->{organismAbbrev};
    my $species = $map->{$organism};
    if (not $species) {
      warn("No production name found for $organism. Do you have a core for this species?\n");
      next;
    }
    
    # Skip if not in the registry to mask
    if (not $to_mask->{$organism}) {
      next;
    }
    
    # Get the reference
    my $ref_organism;
    if ($entry->{is_reference} eq 'TRUE') {
      $ref_organism = $organism;
    } else {
      $ref_organism = $entry->{referenceStrainAbbrev};
    }
    my $ref_species = $map->{$ref_organism};
    if (not $ref_species) {
      warn("No production name found for reference $ref_organism. Do you have a core for this reference?\n");
      next;
    }
    
    # Get repeat library
    my $lib = $ref->{$ref_species};
    if (not $lib) {
      warn("WARNING: No repeat library found for reference $ref_organism ($ref_species). Please check the tab file and update the reference for this species. Skipping this organism.\n");
      next;
    }

    $libraries{$species} = $lib;
  }

  return \%libraries;
}

sub get_libs {
  my ($dir_path) = @_;

  opendir(my $dir, $dir_path);
  my @files = sort grep { not $_ =~ /^\./ } readdir $dir;
  closedir $dir;
  
  # Only keep the filtered version, or unfiltered if there is no filtered version
  my %libs = map { $_ => 1 } @files;
  for my $lib (keys %libs) {
    next if not $lib =~ /\.rm|\.lib|\.filtered/;
    if (exists $libs{$lib . ".filtered"}) {
      delete $libs{$lib};
    }
  }
  
  my %ref_libs;
  for my $lib (keys %libs) {
    my $name = $lib;
    $name =~ s/([^\.]+).*$/$1/;
    $ref_libs{$name} = catfile(rel2abs($dir_path), $lib);
  }
  
  return \%ref_libs;
}

sub get_tab_file {
  my ($tab_path) = @_;

  open my $tab_fh, "<", $tab_path;
  my $head = readline $tab_fh;
  chomp $head;
  my @header = split("\t", $head);
  
  my @entries;
  while (my $line = readline $tab_fh) {
    chomp $line;
    my @values = split "\t", $line;
    my %entry;
    @entry{@header} = @values;
    push @entries, \%entry;
  }
  close $tab_fh;
  
  return \@entries;
}

sub get_species_map {
  my ($registry_path) = @_;

  my $registry = 'Bio::EnsEMBL::Registry';
  my $reg_path = $registry_path;
  $registry->load_all($reg_path);

  my %map;
  my $species = $registry->get_all_species();
  print STDERR scalar(@$species) . " species in registry\n";
  for my $sp (@$species) {
    my $ma = $registry->get_adaptor($sp, "core", "MetaContainer");
    my ($organism_abbrev) = @{ $ma->list_value_by_key('BRC4.organism_abbrev') };
    my ($prod_name) = @{ $ma->list_value_by_key('species.production_name') };
    $map{$organism_abbrev} = $prod_name;
    $ma->dbc->disconnect_if_idle;
  }
  return \%map;
}

sub write_json {
  my ($data, $libs_dir, $out_path) = @_;
  
  open my $json_fh, ">", $out_path;
  print $json_fh JSON->new->allow_nonref->pretty->encode($data);
  close $json_fh;
}

###############################################################################
# Parameters and usage
sub usage {
  my $error = shift;
  my $help = '';
  if ($error) {
    $help = "[ $error ]\n";
  }
  $help .= <<'EOF';
    Get a json parameter file for repeat masking. You need to provide 2 registries:
    * one (ref_registry) with all the species (both prod and new)
    * and one with only the new cores (new_registry)

    --new_registry <path> : Ensembl registry containing all genomes to mask
    --ref_registry <path> : Ensembl registry containing all references
    --tab      <path> : path to the repeat lib file from UPenn
    --libs_dir <path> : libs directory
    --out_json <path> : Json parameter file to create
    
    --help            : show this help message
    --verbose         : show detailed progress
    --debug           : show even more information (for debugging purposes)
EOF
  print STDERR "$help\n";
  exit(1);
}

sub opt_check {
  my %opt = ();
  GetOptions(\%opt,
    "new_registry=s",
    "ref_registry=s",
    "tab=s",
    "libs_dir=s",
    "out_json=s",
    "help",
    "verbose",
    "debug",
  );

  usage("New registry file needed") if not $opt{new_registry};
  usage("Ref registry file needed") if not $opt{ref_registry};
  usage("Tab file needed") if not $opt{tab};
  usage("Libs dir needed") if not $opt{libs_dir};
  usage("Output json needed") if not $opt{out_json};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

