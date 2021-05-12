#!/usr/bin/env perl
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

use strict;
use warnings;
use autodie;

use Carp;
use DateTime;
use Data::Dumper;
use File::Basename;
use File::Path qw(make_path rmtree);
use Storable;

use Bio::EnsEMBL::Utils::CliHelper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

#### get the basic options for connecting to a database server
            #
my $optsd = [
              "input_file:s", 
              "phenotype_dictionary:s", 
              "conditions_dictionary:s",
              "obo_file_destination:s"
            ];

my $opts = $cli_helper->process_args( $optsd, \&pod2usage );

my $phibase_file = $opts->{input_file};
my $phenotype_dictionary = $opts->{phenotype_dictionary};
my $conditions_dictionary = $opts->{conditions_dictionary};
my $obo_file_required;


if ( $opts->{obo_file_destination} ) {
  $obo_file_required = 1;
  print " OBO file is required.\n";
}

print "Phi file: $phibase_file \n";
print "Ontologies phenotypes file: $phenotype_dictionary \n";

##################### List accepted terms ######################################

require $phenotype_dictionary;

print "\n------------- Accepted phenotypes : -------------\n\n";
my %mapped_phenotypes = do $phenotype_dictionary;
my %mapped_phenotypes_values ;

foreach (keys %mapped_phenotypes){
  if (!$mapped_phenotypes_values{$_}){
    $mapped_phenotypes_values{$mapped_phenotypes{$_}} = 1;
  } 
}

my @phenotypes_keys = sort { $a cmp $b } keys %mapped_phenotypes_values; # let's show them nicely!

# foreach my $key ( @phenotypes_keys ) {
#     print "$key \n" ;
# }

print "\n\n------------- Accepted conditions : -------------\n\n";

my %mapped_conditions = do $conditions_dictionary;
my %mapped_conditions_values ;

foreach (keys %mapped_conditions){
  if (!$mapped_conditions_values{$_}){
    $mapped_conditions_values{$mapped_conditions{$_}} = 1;
  } 
}

my @conditions_keys = sort { $a cmp $b } keys %mapped_conditions_values;

# foreach my $key ( @conditions_keys ) {
#     print "$key \n" ;
# }
# print " \n";

########################################################

# open input file
open( my $INPUT_CSV, "<", $phibase_file );

my ($temp_file, $dirs) = fileparse($phibase_file);
my $temp_output_phi_file = "${dirs}Ontologies_temp_files/temp_$temp_file";

# create temp folder
eval { make_path("${dirs}Ontologies_temp_files/" ) };
if ($@) {
  print "Couldn't create ${dirs}Ontologies_temp_files/ : $@";
}

# open and create output temp file
open( my $OUTPUT_CSV, '>' . $temp_output_phi_file);

# my %result_conditions; # list of conditions
# my %result_phenotypes; # list of phenotypes
my %unmapped_phenotypes;
my %unmapped_conditions;


LINE: while (my $line = <$INPUT_CSV> ) {
 
  my @fields = split(/,/, $line);
 

  if (scalar(@fields)>15) {
    print " ********** WARNING: COMMA DETECTED ************ \n" .
          "  A field with an unexpected comma has been found in entry " . $fields[1] . 
          ". This line will be ignored.\n";
    next LINE;
  }

  my $phenotype = lc(remove_spaces($fields[11]));
  my $conditions = lc(remove_spaces($fields[12]));
  my $phi_number = $fields[0];
  my $phi_id = $fields[1];


  my $conditions_field = '';
  my $pubmed_id = $fields[13] || '';
  my $doi = $fields[14] || '';


  if ( $phenotype && $conditions) {
    
    # process phenotypes
    if ($mapped_phenotypes{$phenotype}) { # already exists in the dictionnary
      $phenotype = $mapped_phenotypes{$phenotype};
      # $result_phenotypes{$phenotype}; # add Gene phenotype to results 
    } else {
      if (!$unmapped_phenotypes{$phenotype}) { # if not yet in the list of terms to be solved 
          $unmapped_phenotypes{$phenotype} = $phenotype; # add it 
      }
    } 

    #process conditions
    my @conditions_splitted = split /[;]/, $conditions; # can also try with /[;\/:]/
    my $cond_ctr = 1;
    
    foreach my $cond (@conditions_splitted) {
      $cond = remove_spaces($cond);    
     
      if ($mapped_conditions{$cond}) { # already exists in the dictionnary
        if ($cond_ctr == 1) {
          $conditions_field = $mapped_conditions{$cond};
        } else { # multiple conditions
          $conditions_field =  "$conditions_field; " . $mapped_conditions{$cond};
        }
        #$result_conditions{$cond}; # add to results
      } else {
        if (!$unmapped_conditions{$cond}) { # if not yet in the list of terms to be solved
            $unmapped_conditions{$cond} = $cond; # add it 
        }
      }
      $cond_ctr++;
    }

    my $corrected_line = join(',', $phi_number,
                                $phi_id,
                                $fields[2],
                                $fields[3],
                                $fields[4],
                                $fields[5],
                                $fields[6],
                                $fields[7],
                                $fields[8],
                                $fields[9],
                                $fields[10],
                                $phenotype,
                                $conditions_field,
                                $pubmed_id,
                                $doi);

    print $OUTPUT_CSV $corrected_line;                         
  }
}

close $INPUT_CSV;
close $OUTPUT_CSV;

# print "\n*********** All conditions found in csv file : ************\n";
#  for my $cond ( keys %result_conditions ) {
#    print $cond, "\n";
#  }

# print "\n*********** All phenotypes found in csv file: ************\n";
# for my $phent ( keys %result_phenotypes ) {
#   print $phent, "\n";
# }

print "\n++++++++++++ Unmapped phenotypes : ++++++++++++ \n";
if ( !%unmapped_phenotypes ) {
  print " No unmapped phenotypes found.\n"
} else {
  foreach ( keys %unmapped_phenotypes ) {
    my $new_ph = $unmapped_phenotypes{$_};
    print "\"$new_ph\"\n";
  } 
}

print "\n+++++++++++++ Unmapped conditions : +++++++++++++ \n";
if ( !%unmapped_conditions ) {
  print " No unmapped conditions found.\n"
} else {
  foreach ( keys %unmapped_conditions ) {
    my $new_cond = $unmapped_conditions{$_};
    print "\"$new_cond\"\n";
  } 
}

################ OVERWRITE FILE WITH CORRECT TERMS AND GENERATE OBO FILE ###################

if ( !%unmapped_phenotypes && !%unmapped_conditions ) {
  print "\n>> All ontologies are correctly mapped.\n";
  _overwrite_ontologies();
  if ($obo_file_required) {
    _create_obo_file();
  }
  print "\n";
} else {
  print "\n>> Some unmapped phenotypes or conditions need to be added to its corresponding dictionary.\n"
}

# replace the uncorrect entries in the csv file with the correct terminology 
sub _overwrite_ontologies {
  rename $temp_output_phi_file, $phibase_file;
  rmtree("${dirs}Ontologies_temp_files/"); # remove the temp folder
  print ">> All terms in the original file have been replaced with controlled terms.\n";
};

# create new .obo file with the updated terminology 
sub _create_obo_file { 
  open my $output, ">", $opts->{obo_file_destination} or croak "Could not open " . $opts->{obo_file_destination};
  my $login = getlogin || getpwuid($<) || "User unknown";
  my $date_string = DateTime->now->dmy;
  my $header = "format-version: 5\ndate: $date_string \nsaved-by: $login \n" . 
                "auto-generated-by: construct_phibase_obo.pl\n" .
                "synonymtypedef: systematic_synonym Systematic synonym EXACT\n" .
                "default-namespace: phenotype\nontology: PHI\n\n";

  my $PHENOTYPE_PARENT = 1000000;
  my $CONDITION_PARENT = 2000000;
  my $IDENTIFIER_PARENT = 0;

  #add header
  print $output $header ;
  
  # define parent terms  
  print $output "[Term]\n" .
                "id: PHI:$PHENOTYPE_PARENT\n" .
                "name: phenotype\n" .
                "namespace: phenotype\n" .
                "\n" .
                "[Term]\n" .
                "id: PHI:$CONDITION_PARENT\n" .
                "name: experiment specification\n" .
                "namespace: experiment_specification\n" .
                "\n" . 
                "[Term]\n" .
                "id: PHI:0 \n" .
                "name: phibase identifier\n" .
                "namespace: phibase_identifier\n" .
                "\n"
                ;
  
  # add phenotype terms                    
  my $phi_counter = $PHENOTYPE_PARENT; 
  foreach my $key ( @phenotypes_keys ) {
    if ( $key ne 'phenotype') {
      print $output "[Term]\n" .
                  "id: PHI:" . ++$phi_counter . "\n" .
                  "name: $key\n" .
                  "namespace: phenotype\n" .
                  "is_a: PHI:$PHENOTYPE_PARENT\n\n";
    }
  }
  
  # add condition terms
  $phi_counter = $CONDITION_PARENT; 
  foreach my $key ( @conditions_keys ) {
    if ( $key ne 'experiment specification') {
      print $output "[Term]\n" .
                  "id: PHI:" . ++$phi_counter . "\n" .
                  "name: $key\n" .
                  "namespace: experiment_specification\n" .
                  "is_a: PHI:$CONDITION_PARENT\n\n";
    }
  }
  
  # add identifiers
  $phi_counter = $IDENTIFIER_PARENT; 
  for (1..9999) {
    print $output "[Term]\n" .
                "id: PHI:" . ++$phi_counter . "\n" .
                "name: $phi_counter\n" .
                "namespace: phibase_identifier\n" .
                "is_a: PHI:$IDENTIFIER_PARENT\n\n";
  }
  close $output;
  print ">> A new OBO file has been created in: " . $opts->{obo_file_destination} . "\n";
}

sub remove_spaces {
  my $string = shift;
  $string =~ s/^\s*//;
  $string =~ s/\s*$//;
  return $string;
}
