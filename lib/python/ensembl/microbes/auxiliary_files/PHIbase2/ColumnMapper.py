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

        if 'PHI-base' in db_dict:
            file_rows = {
                    "entry_id": 0,
                    "interactor_A_molecular_id": 1, # uniprot id
                    "interactor_A_sequence": 24,
                    "interactor_A_ensembl_id": 2,
                    "interactor_A_species_taxon_id": 3,
                    "interactor_A_species_strain": 5,
                    "interactor_A_origin_name": 4,  # species name 
                    "interactor_B_molecular_id": 9,   # uniprot id
                    "interactor_B_species_taxon_id": 11,
                    "interactor_B_species_strain": 13,
                    "interactor_B_origin_name": 12,   # species name (equivalent to chemical name when pointing to a synthetic molecule)
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
            ColumnMapper.source_db_name = "PHI-base"
            ColumnMapper.original_curator_db = "PHI-Canto"

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

        if 'PlasticDB' in db_dict:
            file_rows = {
                    "entry_id": 0,
                    "interactor_A_molecular_id": 15,
                    "interactor_A_sequence": 16,
                    "interactor_A_ensembl_id": 14, #uniprot id
                    "interactor_A_species_taxon_id": 2,
                    "interactor_A_species_strain": 28,
                    "interactor_A_origin_name": 1,       #species_name
                    "interactor_B_molecular_id": 4,  #CHEBI or CAS
                    "interactor_B_species_taxon_id": 30,
                    "interactor_B_species_strain": 29,
                    "interactor_B_origin_name": 5,	#chemical name (treated as equivalent of species name)
                    "interactor_B_sequence": 31,
                    "interactor_B_ensembl_id": 32,
                    "litterature_id": 7,
                    "source_db_label": db_dict,
                    }

            ColumnMapper.interactor_A_interactor_type = "protein"
            ColumnMapper.interactor_A_curie_type = "uniprot"
            ColumnMapper.interactor_B_interactor_type = "synthetic"
            ColumnMapper.interactor_B_curie_type = "chebi"

            ColumnMapper.ontology_name = ""
            ColumnMapper.ontology_description = ""
            ColumnMapper.ontology_file = ""

            ColumnMapper.source_db_description = "A database of microorganisms and proteins linked to plastic biodegradation."
            ColumnMapper.source_db_name = "PlasticDB"
            ColumnMapper.original_curator_db = "PlasticDB"

            ColumnMapper.keys_descriptions = {
                    "Enzyme name": "Name of the plastic degradating enzyme",
                    "Experimental evidence": "Experimental evidence",
                    "Plastic used": "Specific plastic used for determining the degradation",
                    "Thermophilic conditions": "Wether thermophilic conditions play a role in the biodegradation",
                    "Isolation environment": "Environment in which the biodegradation was observed",
                    "Isolation location": "Geographical location where the biodegradation was observed",
                    }

            ColumnMapper.keys_rows = {
                    "Enzyme name": 9,
                    "Experimental evidence": 18,
                    "Plastic used": 19,
                    "Thermophilic conditions": 22,
                    "Isolation environment": 24,
                    "Isolation location": 25,
                    }
            ColumnMapper.litterature_source = "DOI"

        if 'HPIDB/(HPIDBcurated)' in db_dict:
            file_rows = {
                    "entry_id": 0,
                    "interactor_A_molecular_id": 1,
                    "interactor_A_sequence": 20,
                    "interactor_A_ensembl_id": 25,
                    "interactor_A_species_taxon_id": 10,
                    "interactor_A_species_strain": 10,
                    "interactor_A_origin_name": 18,    #species name 
                    "interactor_B_molecular_id": 2,
                    "interactor_B_species_taxon_id": 11,
                    "interactor_B_species_strain": 11,
                    "interactor_B_origin_name": 19,   #species name 
                    "interactor_B_sequence": 21,
                    "interactor_B_ensembl_id": 26,
                    "litterature_id": 9,
                    "source_db_label": db_dict,
                    }

            ColumnMapper.interactor_A_interactor_type = "protein"
            ColumnMapper.interactor_A_curie_type = "uniprot"
            ColumnMapper.interactor_B_interactor_type = "protein"
            ColumnMapper.interactor_B_curie_type = "uniprot"

            ColumnMapper.ontology_name = "PSI-MI"
            ColumnMapper.ontology_description = "A structured controlled vocabulary for the annotation of experiments concerned with protein-protein interactions. Developed by the HUPO Proteomics Standards Initiative."
            ColumnMapper.ontology_file = "/nfs/production/flicek/ensembl/microbes/mcarbajo/Phytopath_db/Obo_files/psi-mi.obo"

            ColumnMapper.source_db_description = "A resource that helps annotate, predict and display host-pathogen interactions. https://hpidb.igbb.msstate.edu/index.html"
            ColumnMapper.source_db_name = "HPIDB"
            ColumnMapper.original_curator_db = "HPIDB"

            ColumnMapper.keys_descriptions = {
                    
                    "Experimental evidence": "Detection method or experimental evidence for this interaction",
                    "Interaction type": "Type of interaction",
                    "Confidence": "A method used to derive a numerical or empirical measure of confidence in a particular interaction, or in the identification of the participants in an interaction.",
                    "Interactor_A experimental strain": "Strain in which the interaction was observed (first interactor)",
                    "Interactor_B experimental strain": "Strain in which the interaction was observed (second interactor)",
                    }

            ColumnMapper.keys_rows = {
                    
                    "Experimental evidence": 7,
                    "Interaction type": 12,
                    "Confidence": 14,
                    "Interactor_A experimental strain": 18,
                    "Interactor_B experimental strain": 19,
                    }
            ColumnMapper.litterature_source = "PMID"

        return file_rows
