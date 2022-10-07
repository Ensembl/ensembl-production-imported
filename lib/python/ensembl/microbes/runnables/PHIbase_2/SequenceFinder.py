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
import eHive
import requests
import hashlib

class SequenceFinder(eHive.BaseRunnable):
    """Finds the molecular structure of the interactor"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.warning("Fetch Sequence Finder")
        self.param('failed_job', '')
        phi_id = self.param('PHI_id')
        self.check_param("interactor_A_ensembl_id")
        self.check_param("interactor_B_ensembl_id")
        self.check_param("interactor_A_uniprot_id")
        self.check_param("interactor_B_uniprot_id")

    def run(self):
        self.warning("Sequence finder run")
        self.get_values()

    def get_values(self):
        phi_id = self.param('PHI_id')
        interactor_A_ensembl_gene_stable_id = self.param("interactor_A_ensembl_id")
        interactor_B_ensembl_gene_stable_id = self.param("interactor_B_ensembl_id")
        interactor_A_uniprot_id = self.param("interactor_A_uniprot_id")

        interactor_A_molecular_structure = self.get_molecular_structure(interactor_A_uniprot_id, interactor_A_ensembl_gene_stable_id)
        interactor_B_ensembl_gene_stable_id = self.param('interactor_B_ensembl_id')
        if interactor_B_ensembl_gene_stable_id == "UNDETERMINED":
            interactor_B_molecular_structure = "UNDETERMINED"
        else:
            interactor_B_uniprot_id = self.param("interactor_B_uniprot_id")
            interactor_B_molecular_structure = self.get_molecular_structure(interactor_B_uniprot_id, interactor_B_ensembl_gene_stable_id)
        self.param("interactor_A_molecular_structure",interactor_A_molecular_structure)
        self.param("interactor_B_molecular_structure",interactor_B_molecular_structure)

    def get_molecular_structure(self, uniprot_id, ensembl_gene_id):
        #TO DO: Either redefine with a checksum value of the sequence or remove the sequence completely (we probably don't need it)
        uniprot_seq = 'TO_BE_DEFINED'
        '''
        uniprot_seq = self.get_uniprot_sequence(uniprot_id)
        ensembl_seqs = self.get_ensembl_sequences(ensembl_gene_id)
        phi_id = self.param('PHI_id')
        try:
            if not self.check_equals(uniprot_seq,ensembl_seqs):
                raise (AssertionError)
            else:
                print(f" {phi_id} Sequence match for  uniprot accession {uniprot_id} and ensembl_accession: {ensembl_gene_id}")
        except AssertionError:
            print(f" {phi_id} NO SEQUENCE MATCH for uniprot accession {uniprot_id} and ensembl_accession {ensembl_gene_id}")
        '''
        return uniprot_seq

    def check_equals(self, uniprot_seq, ensembl_seqs):
        for seq in ensembl_seqs:
            if uniprot_seq == seq:
                return True
        return False

    def get_uniprot_sequence(self, uniprot_id):
        uniprot_url = "https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta"
        response = requests.get(uniprot_url)
        uniprot_seq = ''
        for line in response.iter_lines():
            dc_l = line.decode('utf-8')
            if dc_l[0] != ">":
                uniprot_seq = uniprot_seq + str(dc_l)
        return uniprot_seq

    def get_ensembl_sequences(self,ensembl_gene_id):
        ensembl_api_url = "https://rest.ensembl.org/sequence/id/" + ensembl_gene_id + "?type=protein;multiple_sequences=1;content-type=text/x-seqxml%2Bxml"
        response = requests.get(ensembl_api_url, stream=True)
        ensembl_seq_list = []

        for line in response.iter_lines():
            dc_l = line.decode('utf-8').strip().replace("</AAseq>",'')
            if dc_l.startswith("<AAseq>"):
                ensembl_seq_list.append(dc_l.replace("<AAseq>",''))
        return ensembl_seq_list

    def create_digest(input_string):
        digest = hashlib.sha256(str(input_string).encode('utf-8')).hexdigest()
        return digest

    def build_output_hash(self):
        lines_list = []
        entry_line_dict = {
                "interactor_A_molecular_structure": self.param("interactor_A_molecular_structure"),
                "interactor_B_molecular_structure": self.param("interactor_B_molecular_structure"),
                }
        lines_list.append(entry_line_dict)
        return lines_list

    def write_output(self):
        phi_id = self.param('PHI_id')
        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            print(f"{phi_id} written to DBwriter")
            for entry in entries_list:
                self.dataflow(entry, 1)
        else:
            self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))

    def check_param(self, param):
        try:
            self.param_required(param)
        except:
            error_msg = self.param('PHI_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print(error_msg)

