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

# The dictionary is sorted alphabetically by VALUE / KEY

my %mapped_conditions = (
  # PHI-base currently does not treat the ‘experimental conditions’ as controlled terminology and do not filter/check the used terms 
  # but is currently working with the Pombase team to provide a community curation interphase for Pathogen-host interactions. 
  # This work is expected to be presented at the Biocuration Conference in Cambridge, April 2019.
  #
  # This first block use the new ontologies (not completed), the second block uses the old defintions and has a complete dictionary for releases up to 41.
  #
  # NEW ONTOLOGIES (uncomplete. Will work from release 42 on)

  # "altered gene expression / gene regulation"                                       =>  "Altered gene expression / gene regulation",
  # "altered gene expression / gene regulation: overexpression"                       =>  "Altered gene expression / gene regulation: overexpression",
  # "altered gene expression/ gene regulation: overexpression"                        =>  "Altered gene expression / gene regulation: overexpression",
  # "altered gene expression / gene regulation: downregulation"                       =>  "Altered gene expression / gene regulation: downregulation",
  # "altered gene expression / gene regulation: silencing"                            =>  "Altered gene expression / gene regulation: silencing",
  # "altered gene expression/ gene regulation: silencing"                             =>  "Altered gene expression / gene regulation: silencing",
  # "biochemical analysis"                                                            =>  "Biochemical analysis",
  # "biochemical evidence"                                                            =>  "Biochemical analysis",
  # "chemical complementation"                                                        =>  "Chemical complementation",
  # "functional test in host"                                                         =>  "Functional test in host",
  # "functional test in host: dipping"                                                =>  "Functional test in host: dipping",
  # "functional test in host: direct injection"                                       =>  "Functional test in host: direct injection",
  # "functional test in host: transient expression"                                   =>  "Functional test in host: transient expression",
  # "functional test in host: transient expression of recombinant mutant proteins"    =>  "Functional test in host: transient expression of recombinant mutant proteins",
  # "gene complementation"                                                            =>  "Gene complementation",
  # "host gene complementation"                                                       =>  "Gene complementation",
  # "gene deletion: cluster"                                                          =>  "Gene deletion: cluster",
  # "gene deletion: full"                                                             =>  "Gene deletion: full",
  # "gene deletion: partial"                                                          =>  "Gene deletion: partial",
  # "other evidence: host gene disruption"                                            =>  "Gene disruption",
  # "gene disruption"                                                                 =>  "Gene disruption",
  # "other evidence: host gene mutation"                                              =>  "Gene disruption",
  # "gene mutation"                                                                   =>  "Gene mutation",
  # "gene mutation: characterised"                                                    =>  "Gene mutation: characterised",
  # "heterologous expression"                                                         =>  "Heterologous expression",
  # "natural sequence variation"                                                      =>  "Natural sequence variation",
  # "other evidence"                                                                  =>  "Other evidence",
  # "other evidence: antibody inhibition assay"                                       =>  "Other evidence",
  # "other evidence: association genetics"                                            =>  "Other evidence",
  # "other evidence: chemical inhibition"                                             =>  "Other evidence",
  # "other evidence: coinjection"                                                     =>  "Other evidence",
  # "other evidence: co-ip"                                                           =>  "Other evidence",
  # "other evidence: co-ip in vitro"                                                  =>  "Other evidence",
  # "other evidence: enzyme activity inhibition"                                      =>  "Other evidence",
  # "other evidence: gene expression regulation"                                      =>  "Other evidence",
  # "other evidence: insertional mutagenesis screen"                                  =>  "Other evidence",
  # "other evidence: in vitro gst pull-down assay"                                    =>  "Other evidence",
  # "other evidence: mutated transgenic plants"                                       =>  "Other evidence",
  # "other evidence: preventing binding with lipid rafts by depletion of cholesterol" =>  "Other evidence",
  # "other evidence: transgenic plants"                                               =>  "Other evidence",
  # "other evidence: wild-type inoculation"                                           =>  "Other evidence", 
  # "other evidence: yeast split-ubiquitin two-hybrid system"                         =>  "Other evidence", 
  # "other evidence: yeast two-hybrid"                                                =>  "Other evidence",
  # "promoter mutation"                                                               =>  "Promoter mutation",
  # "promotor mutation"                                                               =>  "Promoter mutation",
  # "sequence analysis of sensitive and resistant strains"                            =>  "Sequence analysis of sensitive and resistant strains",
  # "sexual cross, sequencing of resistance conferring allele"                        =>  "Sexual cross, sequencing of resistance conferring allele"

 
#OLD ONTOLOGIES (works for releases up to 41)
  "altered gene expression"                                                          =>  "Altered gene expression",
  "altered gene expression / gene regulation"                                        =>  "Altered gene expression",
  "biochemical analysis"                                                             =>  "Biochemical analysis",
  "biochemical evidence"                                                             =>  "Biochemical analysis",
  "characterised"                                                                    =>  "Characterised gene mutation",
  "characterised gene mutation"                                                      =>  "Characterised gene mutation",
  "gene mutation: characterised"                                                     =>  "Characterised gene mutation",
  "chemical complementation"                                                         =>  "Chemical complementation",
  "cluster deletion"                                                                 =>  "Cluster deletion",
  "direct injection"                                                                 =>  "Direct injection into host organism",
  "functional test in host: direct injection"                                        =>  "Direct injection into host organism",
  "direct injection into host organism"                                              =>  "Direct injection into host organism",
  "other evidence: coinjection"                                                      =>  "Direct injection into host organism",
  "experiment specification"                                                         =>  "experiment specification",
  "full"                                                                             =>  "Full gene deletion",
  "full gene deletion"                                                               =>  "Full gene deletion",
  "gene deletion full"                                                               =>  "Full gene deletion",
  "gene deletion: full"                                                              =>  "Full gene deletion",
  "functional test in host"                                                          =>  "Functional test in host organism",
  "functional test in host organism"                                                 =>  "Functional test in host organism",
  "functional test in host: dipping"                                                 =>  "Functional test in host organism",
  "functional test in host: stable expression of 35s:harxl44 in arabidopsis"         =>  "Functional test in host organism",
  "altered gene expression / gene regulation: silencing: gene complementation"       =>  "Gene complementation",
  "gene complementation"                                                             =>  "Gene complementation",
  "host gene complementation"                                                        =>  "Gene complementation",
  "gene deletion"                                                                    =>  "Gene deletion",
  "gene deletion: partial"                                                           =>  "Gene deletion",
  "other evidence: host gene deletion"                                               =>  "Gene deletion",
  "gene disruption"                                                                  =>  "Gene disruption",
  "other evidence: host gene disruption"                                             =>  "Gene disruption",
  "altered gene expression / gene regulation: downregulation"                        =>  "Gene downregulation",
  "gene downregulation"                                                              =>  "Gene downregulation",
  "downregulation"                                                                   =>  "Gene downregulation",
  "gene mutation"                                                                    =>  "Gene mutation",
  "other evidence: host gene mutation"                                               =>  "Gene mutation",
  "gene silencing"                                                                   =>  "Gene silencing",
  "altered gene expression/ gene regulation: silencing"                              =>  "Gene silencing",
  "altered gene expression / gene regulation: silencing"                             =>  "Gene silencing",
  "silencing"                                                                        =>  "Gene silencing",
  "altered gene expression / gene regulation: overexpression"                        =>  "Gene overexpression",
  "altered gene expression/ gene regulation: overexpression"                         =>  "Gene overexpression",
  "functional test in host: stable overexpression of med19a in arabidopsis"          =>  "Gene overexpression",
  "gene overexpression"                                                              =>  "Gene overexpression",   
  "overexpression"                                                                   =>  "Gene overexpression",
  "heterologous expression"                                                          =>  "Heterologous expression",
  "natural sequence variation"                                                       =>  "Natural sequence variation",
  "partial gene deletion"                                                            =>  "Partial gene deletion",
  "promoter mutation"                                                                =>  "Promoter mutation",
  "sequence analysis"                                                                =>  "Sequence analysis",
  "sequence analysis of sensitive and resistant strains"                             =>  "Sequence analysis of sensitive and resistant strains",
  "sequencing of resistance conferring allele"                                       =>  "Sequencing of resistance conferring allele",
  "sexual cross"                                                                     =>  "Sexual cross",
  "complementation"                                                                  =>  "Other evidence",
  "heterologous complementation"                                                     =>  "Other evidence",
  "other evidence"                                                                   =>  "Other evidence",
  "other evidence: antibody inhibition assay"                                        =>  "Other evidence",
  "other evidence: association genetics"                                             =>  "Other evidence",
  "other evidence: chemical inhibition"                                              =>  "Other evidence",
  "other evidence: co-ip"                                                            =>  "Other evidence",
  "other evidence: co-ip in vitro"                                                   =>  "Other evidence",
  "other evidence: enzyme activity inhibition"                                       =>  "Other evidence",
  "other evidence: gene expression regulation"                                       =>  "Other evidence",
  "other evidence: insertional mutagenesis screen"                                   =>  "Other evidence",
  "other evidence: in vitro gst pull-down assay"                                     =>  "Other evidence",
  "other evidence: mutated transgenic plants"                                        =>  "Other evidence",
  "other evidence: preventing binding with lipid rafts by depletion of cholesterol"  =>  "Other evidence",
  "other evidence: transgenic plants"                                                =>  "Other evidence",
  "other evidence: wild-type inoculation"                                            =>  "Other evidence", 
  "other evidence: yeast split-ubiquitin two-hybrid system"                          =>  "Other evidence", 
  "other evidence: yeast two-hybrid"                                                 =>  "Other evidence",
  "promotor mutation"                                                                =>  "Other evidence",
  "functional test in host: transient expression"                                    =>  "Transient expression in host organism",
  "functional test in host: transient expression of recombinant mutant proteins"     =>  "Transient expression in host organism",
  "transient expression"                                                             =>  "Transient expression in host organism",
  "transient expression in host organism"                                            =>  "Transient expression in host organism",
  "uncharacterised gene mutation"                                                    =>  "Uncharacterised gene mutation",
  "unknown gene expression"                                                          =>  "Unknown gene expression",
);
