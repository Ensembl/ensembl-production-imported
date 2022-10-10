# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


class ColumnMapper():

    ontology_description = ""
    keys_descriptions = {}
    keys_rows = {}
    litterature_source = ""

    def __init__(self, source_db):
        self.__dict__ = self.get_dictionary(source_db)


    def get_dictionary(self, db_dict):
        mapped_dict = {}
        file_rows = {}

        if 'combined_PHIbase' in db_dict:
            file_rows = {
                    "interactor_A_id": 1,
                    "interactor_A_sequence": 24,
                    "interactor_A_ensembl_id": 2,
                    "interactor_A_species_taxon_id": 3,
                    "interactor_A_species_strain": 5,
                    "interactor_A_name": 4,
                    "interactor_B_id": 9,
                    "interactor_B_species_taxon_id": 11,
                    "interactor_B_species_strain": 13,
                    "interactor_B_name": 12,
                    "interactor_B_sequence": 25,
                    "interactor_B_ensembl_id": 10,
                    "litterature_id": 22,
                    "source_db_label": db_dict,
                    }

            ColumnMapper.interactor_A_interactor_type = "protein"
            ColumnMapper.interactor_A_curie_type = "uniprot"
            ColumnMapper.interactor_B_interactor_type = "protein"
            ColumnMapper.interactor_B_curie_type = "uniprot"

            ColumnMapper.ontology_description = "Expertly curated molecular and biological information on genes proven to affect the outcome of pathogen-host interactions"

            ColumnMapper.ontology_file = "/hps/software/users/ensembl/repositories/microbes/mcarbajo/checkout/ensembl-production-imported/lib/python/ensembl/microbes/auxiliary_files/PHIbase2/PHI_DUMMY.obo"
                    
            ColumnMapper.source_db_description = "Pathogen-Host Interactions Database that catalogues experimentally verified pathogenicity."
            
            ColumnMapper.keys_descriptions = {
                    "Interaction phenotype": "Interaction phenotype/ disease outcome",
                    "Disease name": "Name of disease",
                    "Pathogen protein modification:": "Protein modification in the pathogen interactor",
                    "Host protein modification": "Protein modification in the host interactor",
                    "Experimental evidence": "Experimental evidence",
                    "Interaction type": "Type of interaction",
                    "PHI-base high level term": "High level phenotype category used by PHI-base to characterise the mutant",
                    "Pathogen experimental strain": "Pathogen strain in which the interaction was observed",
                    "Host experimental strain": "Host strain in which the interaction was observed", 
                    }

            ColumnMapper.keys_rows = {
                    "Interaction phenotype": 17,
                    "Disease name": 20,
                    "Pathogen protein modification": 8,
                    "Host protein modification": 16,
                    "Experimental evidence": 20,
                    "Interaction type": 21,
                    "PHI-base high level term": 23,
                    "Pathogen experimental strain":6,
                    "Host experimental strain":14,
                    }
            ColumnMapper.litterature_source = "PMID"
        return file_rows
