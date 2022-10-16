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
                    "interactor_A_molecular_id": 1,
                    "interactor_A_sequence": 24,
                    "interactor_A_ensembl_id": 2,
                    "interactor_A_species_taxon_id": 3,
                    "interactor_A_species_strain": 5,
                    "interactor_A_name": 4,  #species name 
                    "interactor_B_molecular_id": 9,
                    "interactor_B_species_taxon_id": 11,
                    "interactor_B_species_strain": 13,
                    "interactor_B_name": 12,   #species name or chemical name if it points to a synthetic molecule
                    "interactor_B_sequence": 25,
                    "interactor_B_ensembl_id": 10,
                    "litterature_id": 22,
                    "source_db_label": db_dict,
                    }

            ColumnMapper.interactor_A_interactor_type = "protein"
            ColumnMapper.interactor_A_curie_type = "uniprot"
            ColumnMapper.interactor_B_interactor_type = "protein"
            ColumnMapper.interactor_B_curie_type = "uniprot"

            ColumnMapper.ontology_name = "PHIPO"
            ColumnMapper.ontology_description = "Pathogen Host Interactions Phenotype Ontology. Ontology of species-neutral phenotypes observed in pathogen-host interactions."

            ColumnMapper.ontology_file = "/nfs/production/flicek/ensembl/microbes/mcarbajo/Phytopath_db/Obo_files/phipo-simple.obo"
                    
            ColumnMapper.source_db_description = "Pathogen-Host Interactions Database that catalogues experimentally verified pathogenicity."
            
            ColumnMapper.keys_descriptions = {
                    "Interaction phenotype": "Interaction phenotype/ disease outcome",
                    "Disease name": "Name of disease",
                    "Pathogen protein modification": "Protein modification in the pathogen interactor",
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

        if 'plasticDB' in db_dict:
            file_rows = {
                    "interactor_A_id": 14,
                    "interactor_A_sequence": 15,
                    "interactor_A_ensembl_id": 13,
                    "interactor_A_species_taxon_id": 1,
                    "interactor_A_species_strain": 27,
                    "interactor_A_name": 0,
                    "interactor_B_id": 3,
                    "interactor_B_species_taxon_id": 29,
                    "interactor_B_species_strain": 28,
                    "interactor_B_name": 4,			#species name or chemical name if it points to a synthetic molecule
                    "interactor_B_sequence": 30,
                    "interactor_B_ensembl_id": 31,
                    "litterature_id": 6,
                    "source_db_label": db_dict,
                    }

            ColumnMapper.interactor_A_interactor_type = "protein"
            ColumnMapper.interactor_A_curie_type = "uniprot"
            ColumnMapper.interactor_B_interactor_type = "synthetic"
            ColumnMapper.interactor_B_curie_type = "chebi"

            ColumnMapper.ontology_description = ""

            ColumnMapper.ontology_file = ""

            ColumnMapper.source_db_description = "A database of microorganisms and proteins linked to plastic biodegradation."

            ColumnMapper.keys_descriptions = {
                    "Enzyme name": "Name of the plastic degradating enzyme",
                    "Experimental evidence": "Experimental evidence",
                    "Plastic used": "Specific plastic used for determining the degradation",
                    "Thermophilic conditions": "Wether thermophilic conditions play a role in the biodegradation",
                    "Isolation environment": "Environment in which the biodegradation was observed",
                    "Isolation location": "Geographical location where the biodegradation was observed",
                    }

            ColumnMapper.keys_rows = {
                    "Enzyme name": 8,
                    "Experimental evidence": 17,
                    "Plastic used": 18,
                    "Thermophilic conditions": 21,
                    "Isolation environment": 23,
                    "Isolation location": 24,
                    }
            ColumnMapper.litterature_source = "DOI"

        return file_rows
