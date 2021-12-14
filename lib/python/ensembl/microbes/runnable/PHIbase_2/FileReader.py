
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

# standaloneJob.pl eHive.examples.TestRunnable -language python3

import os
import subprocess
import unittest
import csv
import codecs
import eHive

class FileReader(eHive.BaseRunnable):
    """csv parser to fan 1 job per column"""

    def param_defaults(self):
        return {
            'inputfile' : '#inputfile#',
            'registry'  : '#registry#',
            'branch_to_flow_on_fail' :  -3,
        }

    def fetch_input(self):
        print("inputfile is", self.param_required('inputfile'))
        print("registry", self.param_required('registry'))

    def run(self):
        self.param('entries_list', self.read_lines())

    def write_output(self):
        entries_list = self.param('entries_list')
        for entry in entries_list:
            self.dataflow(entry, 2)

    def read_lines(self):
        self.warning("read_lines running")
        int_db_url, ncbi_tax_url, meta_db_url = self.read_registry()
        self.param("interactions_db_url",int_db_url)
        with open(self.param('inputfile'), newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            next(reader)
            lines_list = []
            for row in reader:
                entry_line_dict = {
                    "branch_to_flow_on_fail" : self.param('branch_to_flow_on_fail'),
                    "PHI_id": row[0],
                    "patho_uniprot_id": row[2],
                    "patho_sequence": row[5],
                    "patho_interactor_name": row[8],
                    "patho_protein_modification": self.limit_string_length(row[10]),
                    "patho_other_names": row[4],
                    "patho_species_taxon_id": row[15],
                    "patho_species_strain": row[17],
                    "host_uniprot_id": self.get_uniprot_id(row[47]),
                    "host_interactor_name": row[46],
                    "host_other_names": row[47],
                    "host_species_taxon_id": row[21],
                    "host_protein_modification": self.limit_string_length(row[48]),
                    "disease_name": row[19],
                    "interaction_phenotype": self.limit_string_length(row[50]),
                    "litterature_id": row[56],
                    "litterature_source": row[57],
                    "doi": row [58],
                    "pathogen_mutant_phenotype": self.limit_string_length(row[32]),
                    "experimental_evidence": self.limit_string_length(row[52]),
                    "transient_assay_exp_ev": self.limit_string_length(row[53]),
                    "host_response_to_pathogen": self.limit_string_length(row[51]),
                    "interactions_db_url": int_db_url,
                    "ncbi_taxonomy_url": ncbi_tax_url,
                    "meta_ensembl_url": meta_db_url,
                    "source_db_label": self.get_db_label(),
        }   
                lines_list.append(entry_line_dict)
            return lines_list

    def get_uniprot_id(self,accessions):
        result = None
        reported = False
        accessions_list = accessions.split(';')
        for ac in accessions_list:
            if "Uniprot: " in ac:
                result = ac.replace('Uniprot: ','')
                reported = True
        try:
            if not reported:
                raise (AssertionError)
        except AssertionError:
            print("Mmmm..., we need to find out the uniprot accession for this interactor: " + accessions)
        return result

    def get_db_label(self):
	#TODO: implement a discriminating method to determine provenance of file(PHI-base, IntAct, HPIDD...)
        return 'PHI-base'
    
    def limit_string_length(self, string):
        value = "%.255s" % string
        print (f"limited value={value}")
        return value

    def read_registry(self):
        with open(self.param('registry'), newline='') as reg_file:
            url_reader = csv.reader(reg_file, delimiter='\t')
            for url in url_reader:
                if url[0] == 'interactions_db_url':
                    int_db_url=url[1]
                elif url[0] == 'ncbi_tax_url':
                    ncbi_tax_url=url[1]
                elif url[0] == 'meta_db_url':
                    meta_db_url=url[1]

            return int_db_url, ncbi_tax_url, meta_db_url

