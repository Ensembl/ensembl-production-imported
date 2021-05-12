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

my %mapped_phenotypes = (
  "phenotype"                                   =>  "phenotype",
  "chemistry target: phenotype unknown"         =>  "Chemistry target: phenotype unknown",
  "chemistry target: resistance to chemical"    =>  "Chemistry target: resistance to chemical",
  "chemistry target: sensitivity to chemical"   =>  "Chemistry target: sensitivity to chemical",
  "effector (plant avirulence determinant)"     =>  "Effector (plant avirulence determinant)",
  "enhanced antagonism"                         =>  "Enhanced antagonism",
  "increased virulence"                         =>  "Increased virulence",
  "increased virulence (hypervirulence)"        =>  "Increased virulence (hypervirulence)",
  "lethal"                                      =>  "Lethal",
  "loss of pathogenicity"                       =>  "Loss of pathogenicity",
  "mixed outcome"                               =>  "Mixed outcome",
  "reduced virulence"                           =>  "Reduced virulence",
  "unaffected pathogenicity"                    =>  "Unaffected pathogenicity", 
  "wild-type mutualism"                         =>  "Wild-type mutualism",
);
