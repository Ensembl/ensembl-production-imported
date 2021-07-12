
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
    """csv parser to fan 1 job per column0"""

    def param_defaults(self):
        return {
            'inputfile' : '#inputfile#',
            'interactions_db_url' : '#interactions_db_url#',
        }

    def fetch_input(self):
        self.warning("Fetch input!")
        print("inputfile is", self.param_required('inputfile'))
        print("interactions_db_url is ", self.param_required('interactions_db_url'))
    
    def run(self):
        self.warning("FileReader run")
        self.param('entries_list', self.read_lines())

    def write_output(self):
        self.warning("Write to the world !")
        entries_list = self.param('entries_list')
        for entry in entries_list:
            self.dataflow(entry, 2)

    def read_lines(self):
        self.warning("read_lines running")
        int_db_url = self.param_required('interactions_db_url')
        with open(self.param('inputfile'), newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            next(spamreader)
            lines_list = []
            for row in spamreader:
                entry_line_dict = {
                    "PHI_id": row[0],
                    "patho_uniprot_id": row[2],
                    "patho_sequence": row[5],
                    "patho_gene_name": row[8],
                    "patho_specie_taxon_id": row[15],
                    "patho_species_name": row[16],
                    "patho_species_strain": row[17],
                    "host_uniprot_id": row[47],
                    "host_gene_name": row[46],
                    "host_species_taxon_id": row[21],
                    "host_species_name": row[22],
                    "litterature_id": row[56],
                    "litterature_source": row[57],
                    "doi": row [58],
                    "interaction_phenotype": row[50],
                    "pathogen_mutant_phenotype": row[32],
                    "experimental_evidence": row[52],
                    "transient_assay_exp_ev": row[53],
                    "interactions_db_url": int_db_url,
                }   
                lines_list.append(entry_line_dict)
            return lines_list

