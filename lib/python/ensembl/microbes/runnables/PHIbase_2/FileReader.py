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

import os
import subprocess
import unittest
import csv
import codecs
import eHive
import ensembl.microbes.auxiliary_files.PHIbase2.ColumnMapper as col_map

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
        int_db_url, ncbi_tax_url, meta_db_url, vertebrate_url, non_vertebrate_url, bacteria_url  = self.read_registry()
        self.param("interactions_db_url",int_db_url)
        
        source_db = self.param('source_db')
        cm = col_map.ColumnMapper(source_db)
        print (f"COLUMN Mapper cm.source_db_label {source_db}")
        with open(self.param('inputfile'), newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            next(reader)
            lines_list = []
            for row in reader:
                try:
                    entry_line_dict = {
                        "branch_to_flow_on_fail" : self.param('branch_to_flow_on_fail'),
                        "entry_id": row[cm.entry_id],
                        "interactor_A_molecular_id": row[cm.interactor_A_molecular_id],#either uniprot or chebi
                        "interactor_A_interactor_type": cm.interactor_A_interactor_type,
                        "interactor_A_curie_type": cm.interactor_A_curie_type,
                        "interactor_A_sequence": row[cm.interactor_A_sequence],
                        "interactor_A_ensembl_id": row[cm.interactor_A_ensembl_id],
                        "interactor_A_species_taxon_id": row[cm.interactor_A_species_taxon_id],
                        "interactor_A_species_strain": row[cm.interactor_A_species_strain],
                        "interactor_A_origin_name": row[cm.interactor_A_origin_name], #either species name or chebi
                        "interactor_B_molecular_id": row[cm.interactor_B_molecular_id],#either uniprot or chebi
                        "interactor_B_interactor_type": cm.interactor_B_interactor_type,
                        "interactor_B_curie_type": cm.interactor_B_curie_type,
                        "interactor_B_sequence": row[cm.interactor_B_sequence],
                        "interactor_B_ensembl_id": row[cm.interactor_B_ensembl_id],
                        "interactor_B_species_taxon_id": row[cm.interactor_B_species_taxon_id],
                        "interactor_B_species_strain": row[cm.interactor_B_species_strain],
                        "interactor_B_origin_name": row[cm.interactor_B_origin_name], #either species name or chebi
                        "litterature_id": row[cm.litterature_id],
                        "litterature_source": cm.litterature_source,
                        "doi": row [cm.litterature_id],
                        "interactions_db_url": int_db_url,
                        "ncbi_taxonomy_url": ncbi_tax_url,
                        "meta_ensembl_url": meta_db_url,
                        "vertebrate_url":vertebrate_url,
                        "non_vertebrate_url": non_vertebrate_url,
                        "bacteria_url":bacteria_url,
                        "source_db_label": cm.source_db_label,
                        "source_db_description": cm.source_db_description,
                        "obo_file": cm.ontology_file,
                        }
                    keys_rows_dict = col_map.ColumnMapper.keys_rows
                    for key in keys_rows_dict:
                        row_number = keys_rows_dict[key]                  
                        row_entry = row[row_number]
                        if row_entry != '':
                            entry_line_dict[key] = self.limit_string_length(row_entry)
                    lines_list.append(entry_line_dict)
                except Exception as e:
                    print(e)
                    pass
        return lines_list


    def get_db_label(self):
        source_db = self.param('source_db')
        return source_db
    
    def limit_string_length(self, string):
        value = "%.255s" % string
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
                elif url[0] == 'vertebrate_url': 
                    vertebrate_url=url[1]
                elif url[0] == 'non_vertebrate_url': 
                    non_vertebrate_url=url[1] 
                elif url[0] == 'bacteria_url': 
                    bacteria_url=url[1] 
                    
        return int_db_url, ncbi_tax_url, meta_db_url,vertebrate_url, non_vertebrate_url, bacteria_url

