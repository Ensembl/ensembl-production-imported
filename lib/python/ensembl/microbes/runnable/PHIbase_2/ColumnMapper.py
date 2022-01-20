#Ensembl
#Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2021] EMBL-European Bioinformatics Institute
#
#This product includes software developed at:
#
#EMBL-European Bioinformatics Institute
#Wellcome Trust Sanger Institute

class ColumnMapper():

    ontology_description = ""
    keys_descriptions = {}
    keys_rows = {}

    def __init__(self, source_db):
        self.__dict__ = self.get_dictionary(source_db)


    def get_dictionary(self, d):
        mapped_dict = {}
        if (d == 'PHI-base'):
            file_rows = {
                    "patho_uniprot_id": 2,
                    "patho_sequence": 5,
                    "patho_interactor_name": 8,
                    "patho_other_names": 4,
                    "patho_species_taxon_id": 15,
                    "patho_species_strain": 17,
                    "host_uniprot_id": 47,
                    "host_interactor_name": 46,
                    "host_other_names": 47,
                    "host_species_taxon_id": 21,
                    "interaction_phenotype": 50,
                    "litterature_id": 56,
                    "litterature_source": 57,
                    "doi": 58,
                    "pathogen_mutant_phenotype": 32,
                    "source_db_label": d,
                    }

            ColumnMapper.ontology_description = "Expertly curated molecular and biological information on genes proven to affect the outcome of pathogen-host interactions"

            ColumnMapper.keys_descriptions = {
                    "interaction_phenotype": "interaction phenotype/ disease outcome",
                    "disease_name": "name of disease",
                    "patho_protein_modification": "protein modification in the pathogen interactor",
                    "host_protein_modification": "protein modification in the host interactor",
                    "experimental_evidence": "experimental evidence",
                    "transient_assay_exp_ev": "transient experimental evidence",
                    "host_response_to_pathogen": "host response to the pathogen",
                    }

            ColumnMapper.keys_rows = {
                    "interaction_phenotype": 50,
                    "disease_name": 19,
                    "patho_protein_modification": 10,
                    "host_protein_modification": 48,
                    "experimental_evidence": 52,
                    "transient_assay_exp_ev": 53,
                    "host_response_to_pathogen": 51,
                    }

        return file_rows
