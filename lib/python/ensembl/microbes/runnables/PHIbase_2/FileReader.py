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
        int_db_url, ncbi_tax_url, meta_db_url = self.read_registry()
        self.param("interactions_db_url",int_db_url)
        
        source_db = self.param('source_db')
        cm = col_map.ColumnMapper(source_db)
        print (f"COLUMN Mapper cm.source_db_label {cm.source_db_label}")
        with open(self.param('inputfile'), newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            next(reader)
            lines_list = []
            for row in reader:
                entry_line_dict = {
                    "branch_to_flow_on_fail" : self.param('branch_to_flow_on_fail'),
                    "PHI_id": row[0],
                    "patho_uniprot_id": row[cm.patho_uniprot_id],
                    "patho_sequence": row[cm.patho_sequence],
                    "patho_interactor_name": row[cm.patho_interactor_name],
                    "patho_other_names": row[cm.patho_other_names],
                    "patho_species_taxon_id": row[cm.patho_species_taxon_id],
                    "patho_species_strain": row[cm.patho_species_strain],
                    "host_uniprot_id": self.get_uniprot_id(row[cm.host_uniprot_id]),
                    "host_interactor_name": row[cm.host_interactor_name],
                    "host_other_names": row[cm.host_other_names],
                    "host_species_taxon_id": row[cm.host_species_taxon_id],
                    "litterature_id": row[cm.litterature_id],
                    "litterature_source": row[cm.litterature_source],
                    "doi": row [cm.doi],
                    "interactions_db_url": int_db_url,
                    "ncbi_taxonomy_url": ncbi_tax_url,
                    "meta_ensembl_url": meta_db_url,
                    "source_db_label": cm.source_db_label,
                    }
                keys_rows_dict = col_map.ColumnMapper.keys_rows
                for key in keys_rows_dict:
                    row_number = keys_rows_dict[key]                  
                    row_entry = row[row_number]
                    if row_entry != '':
                        entry_line_dict[key] = self.limit_string_length(row_entry)
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
            print("Uniprot id not found. Thi interactor needs an Uniprot accession: " + accessions)
        return result

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

            return int_db_url, ncbi_tax_url, meta_db_url

