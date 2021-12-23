# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2021] EMBL-European Bioinformatics Institute
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

class ColumnMapper():

    def __init__(self, d):
        self.__dict__ = self.get_dictionary(d)


    def get_dictionary(self, d):
        mapped_dict = {}
        if (d == 'PHI-base'):
            mapped_dict = {
                "patho_uniprot_id": 2,
                "patho_sequence": 5,
                "patho_interactor_name": 8,
                "patho_protein_modification": 10,
                "patho_other_names": 4,
                "patho_species_taxon_id": 15,
                "patho_species_strain": 17,
                "host_uniprot_id": 47,
                "host_interactor_name": 46,
                "host_other_names": 47,
                "host_species_taxon_id": 21,
                "host_protein_modification": 48,
                "disease_name": 19,
                "interaction_phenotype": 50,
                "litterature_id": 56,
                "litterature_source": 57,
                "doi": 58,
                "pathogen_mutant_phenotype": 32,
                "experimental_evidence": 52,
                "transient_assay_exp_ev": 53,
                "host_response_to_pathogen": 51,
                "source_db_label": d,
                }
        return mapped_dict