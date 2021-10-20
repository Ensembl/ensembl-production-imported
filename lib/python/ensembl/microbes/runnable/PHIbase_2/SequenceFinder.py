
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
import sqlalchemy as db
import sqlalchemy_utils as db_utils
import pymysql
import eHive
import datetime
import re
import ensembl.microbes.runnable.PHIbase_2.core_DB_models as core_db_models
import ensembl.microbes.runnable.PHIbase_2.interaction_DB_models as interaction_db_models
import requests
from xml.etree import ElementTree
import time
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

pymysql.install_as_MySQLdb()

class SequenceFinder(eHive.BaseRunnable):
    """Finds the molecular structure of the interactor"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.warning("Fetch Sequence Finder")
        phi_id = self.param_required('PHI_id')
        print(f'phi_id--{phi_id}')
        core_db_url = self.param_required('core_db_url')
        (c_user,c_pwd,c_host,c_port,c_db) = re.compile(jdbc_pattern).findall(core_db_url)[0]
        self.param('core_user',c_user)
        self.param('core_pwd',c_pwd)
        self.param('core_host',c_host)
        self.param('core_port',int(c_port))
        self.param('core_db',c_db)


    def run(self):
        self.warning("Sequence finder run")
        self.get_values()

    def get_values(self):

        phi_id = self.param_required('PHI_id')
        patho_ensembl_gene_stable_id = self.param_required("patho_ensembl_gene_stable_id")
        host_ensembl_gene_stable_id= self.param_required("host_ensembl_gene_stable_id")
        patho_uniprot_id = self.param_required("patho_uniprot_id")
        host_uniprot_id = self.param_required("host_uniprot_id")

        patho_ensembl_gene_stable_id = self.param_required("patho_ensembl_gene_stable_id")
        host_ensembl_gene_stable_id = self.param_required("host_ensembl_gene_stable_id")
        
        patho_molecular_structure = self.get_molecular_structure(self.param("patho_uniprot_id"), patho_ensembl_gene_stable_id)
        host_molecular_structure = self.get_molecular_structure(self.param("host_uniprot_id"), host_ensembl_gene_stable_id)
        self.param("patho_molecular_structure",patho_molecular_structure)
        self.param("host_molecular_structure",host_molecular_structure)

    def get_molecular_structure(self, uniprot_id, ensembl_gene_id):
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

    def build_output_hash(self):
        lines_list = []
        entry_line_dict = {
                "patho_molecular_structure": self.param("patho_molecular_structure"),
                "host_molecular_structure": self.param("host_molecular_structure"),
                }
        lines_list.append(entry_line_dict)
        return lines_list

    def write_output(self):
        entries_list = self.build_output_hash()
        for entry in entries_list:
            self.dataflow(entry, 1)
