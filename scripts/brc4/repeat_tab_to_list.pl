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
use Capture::Tiny qw(capture);
use Log::Log4perl qw( :easy );
Log::Log4perl->easy_init($WARN);
my $logger = get_logger();
use File::Spec::Functions qw(rel2abs catfile);

use Bio::EnsEMBL::Registry;

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $ref_libs = get_libs($opt{libs_dir});
my $tab = get_references($opt{tab});
my $map = get_ref_map($opt{ref_registry}, $opt{verbose});

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($opt{registry}, 1);

my $sps = $registry->get_all_species();
print STDERR scalar(@$sps) . " cores\n";
die "No core found in the registry" if @$sps == 0;

my %new_map;
for my $species (sort @$sps) {
  my $meta = $registry->get_adaptor($species, 'core', 'MetaContainer');
  if (not $meta) {
    warn "Can't find meta for the core of $species";
  }

  my ($org) = @{ $meta->list_value_by_key("brc4.organism_abbrev") };
  my ($db_gca) = @{ $meta->list_value_by_key("assembly.accession") };

  # We only need to consider the cores for the reference species
  if (exists $tab->{$org}) {
    my $ref_org = $tab->{$org};

    if (exists $map->{$ref_org}) {
      my $ref_name = $map->{$ref_org};
      $new_map{$ref_name} = $ref_org;
    } else {
      warn("No reference core genome found for $ref_org (reference for $org). Using itself.\n");
      $new_map{$species} = $org;
    }
  } else {
    warn("The organism $org is not in the tab file\n");
  }
  $meta->dbc->disconnect_if_idle();
}

# Output the list, filtered
my $outfh;
if ($opt{output_to_model}) {
  open $outfh, ">", $opt{output_to_model};
} else {
  $outfh = *STDOUT;
}
my $nlibs = 0;
for my $map_ref (sort keys %new_map) {
  next if exists $ref_libs->{$map_ref};
  print $outfh "$map_ref\t$new_map{$map_ref}\n";
  $nlibs++;
}
close($outfh) if $opt{output_to_model};
print STDERR "$nlibs new libraries to build\n";

###############################################################################
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

sub get_ref_map {
  my ($reg_path, $verbose) = @_;

  my %map;
  my ($stdout, $stderr, $status) = capture {
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_all($reg_path, 1);

    my $sps = $registry->get_all_species();
    for my $species (sort @$sps) {
      my $meta = $registry->get_adaptor($species, 'core', 'MetaContainer');
      if (not $meta) {
        warn "Can't find meta for the core of $species";
      }

      my ($org) = @{ $meta->list_value_by_key("brc4.organism_abbrev") };
      my ($prod_name) = @{ $meta->list_value_by_key("species.production_name") };
      $map{$org} = $prod_name;
      $meta->dbc->disconnect_if_idle();
    }
  };
  if ($stderr) {
    if ($verbose) {
      warn "Registry loading full log: $stderr\n";
    } else {
      my $nlines = scalar(split "\n", $stderr);
      warn "The Registry loading sent errors ($nlines lines), usually just release differences (check with -v)\n";
    }
  }
  return \%map;
}

sub get_references {
  my ($tab_path) = @_;

  my %data;

  # First, extract the data from the tab file
  open my $tab_fh, "<", $tab_path;
  my @head = split /\t/, readline $tab_fh;
  while (my $line = readline $tab_fh) {
    chomp $line;
    my @values = split /\t/, $line;
    map { $_ =~ s/^ +//; $_ =~ s/ +$// } @values;
    my %sp_data;
    @sp_data{@head} = @values;

    #    next if not $sp_data{"INSDC Accession"};
    if (not ($sp_data{"organismAbbrev"})) {
      warn("Missing organism_abbrev in: $line");
    } else {
      $data{$sp_data{organismAbbrev}} = \%sp_data;
    }
  }
  close $tab_fh;

  # Check that the reference exists (and has an accession) and store it with its accession
  my %ref_data;
  for my $sp (keys %data) {
    my $ref_sp = $data{$sp}->{referenceStrainAbbrev};
    if (not $ref_sp) {
      warn("No reference species given for $sp");
    }
    elsif (not exists $data{$ref_sp}) {
      warn("The reference species ($ref_sp) for species $sp is not in the tab\n");
    }
    $ref_data{$sp} = $ref_sp;
  }

  return \%ref_data;
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
    Get a list of production_names from a repeat tab from BRC4

    --tab <path>      : repeat tab file
    --registry <path> : Ensembl registry for the new genomes
    --ref_registry <path> : Ensembl registry for the reference genomes
    --libs_dir <path> : Path to existing repeat libraries
    --output_to_model <path> : List of species to provide to the modeler

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
    "tab=s",
    "registry=s",
    "ref_registry=s",
    "libs_dir=s",
    "output_to_model=s",
    "help",
    "verbose",
    "debug",
  );

  usage("tab file needed") if not $opt{tab};
  usage("registry needed") if not $opt{registry};
  usage("ref_registry needed") if not $opt{ref_registry};
  usage("libs_dir needed") if not $opt{libs_dir};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

