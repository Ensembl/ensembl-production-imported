
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
        }

    def fetch_input(self):
        self.warning("Fetch input!")
        print("inputfile is", self.param_required('inputfile'))
        print("registry", self.param_required('registry'))

    def run(self):
        self.warning("FileReader run")
        self.param('entries_list', self.read_lines())

    def write_output(self):
        entries_list = self.param('entries_list')
        for entry in entries_list:
            self.dataflow(entry, 2)

    def read_lines(self):
        self.warning("read_lines running")
        int_db_url, ncbi_tax_url, meta_db_url = self.read_registry()
        print("int_db_url:", int_db_url)
        print("ncbi_tax_url:", ncbi_tax_url)
        print("meta_db_url:", meta_db_url)
        with open(self.param('inputfile'), newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            next(reader)
            lines_list = []
            for row in reader:
                entry_line_dict = {
                    "PHI_id": row[0],
                    "patho_uniprot_id": row[2],
                    "patho_sequence": row[5],
                    "patho_gene_name": row[4],
                    "patho_species_taxon_id": row[15],
                    "patho_species_strain": row[17],
                    "host_uniprot_id": row[47],
                    "host_gene_name": row[47],
                    "host_species_taxon_id": row[21],
                    "litterature_id": row[56],
                    "litterature_source": row[57],
                    "doi": row [58],
                    "interaction_phenotype": row[50],
                    "pathogen_mutant_phenotype": row[32],
                    "experimental_evidence": row[52],
                    "transient_assay_exp_ev": row[53],
                    "interactions_db_url": int_db_url,
                    "ncbi_taxonomy_url": ncbi_tax_url,
                    "meta_ensembl_url": meta_db_url,
                }   
                lines_list.append(entry_line_dict)
            return lines_list

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

