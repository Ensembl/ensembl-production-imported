
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
import sqlalchemy as db
import sqlalchemy_utils as db_utils
import pymysql
import eHive
import datetime
import re
import models as core_db_models
import ensembl.microbes.auxiliary_files.PHIbase2.interaction_DB_models as interaction_db_models
import requests
from xml.etree import ElementTree
import time
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

pymysql.install_as_MySQLdb()

class EnsemblCoreReader(eHive.BaseRunnable):
    """Reads from CoreDB to gather fields to attach interactor to the ensembl_gene (mostly ensembl_gene_stable_id)"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.param('failed_job', '')
        source_db = self.param("source_db_label")
        
        if "PHI-base" in source_db:
            phi_id = self.param_required('PHI_id')

        self.check_param('patho_division')
        self.check_param('host_division')
        self.check_param('patho_species_name')
        self.check_param('patho_species_taxon_id')
        self.check_param('host_species_name')
        self.check_param('host_species_taxon_id')
        self.check_param('patho_core_dbname')
        self.check_param('host_core_dbname')
    
    def run(self):
        self.warning("EnsemblCoreReader run")
        self.get_values()

    def get_values(self):     
        patho_species_taxon_id = int(self.param_required('patho_species_taxon_id'))
        patho_division = self.param_required('patho_division')
        host_species_taxon_id = int(self.param_required('host_species_taxon_id'))
        host_division = self.param_required('host_division')
        
       
        try:
            patho_ensembl_gene_stable_id = self.param('patho_ensembl_id')
        except Exception :
            patho_ensembl_gene_stable_id = self.get_ensembl_id(patho_species_taxon_id, self.param("patho_uniprot_id"), self.param("patho_production_name"))

        try:    
            host_ensembl_gene_stable_id = self.param('host_ensembl_id')
        except Exception:
            try:
                host_uniprot_id = self.param("host_uniprot_id")
                host_ensembl_gene_stable_id = self.get_ensembl_id(host_species_taxon_id, host_uniprot_id, self.param("host_production_name"))
            except Exception:
                host_ensembl_gene_stable_id = "UNDETERMINED"

        self.param("patho_ensembl_gene_stable_id",patho_ensembl_gene_stable_id)
        self.param("host_ensembl_gene_stable_id",host_ensembl_gene_stable_id)
    

    def get_ensembl_id(self, tax_id, uniprot_id, species_production_name):
        print("... species_prod_name" + species_production_name)
        url = "https://rest.ensembl.org/xrefs/symbol/" + species_production_name + "/" + uniprot_id + "?external_db=UNIPROT;content-type=application/json"
        print("**-** url" + url)
        response = requests.get(url)
        print("** response" + str(response))


    def get_ensembl_gene_value(self, session, stable_id, species_id):
        try:
            ensembl_gene_value = session.query(interaction_db_models.EnsemblGene).filter_by(ensembl_stable_id=stable_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for species: {species_id} - gene stable_id {stable_id}")
        except NoResultFound:
            ensembl_gene_value = interaction_db_models.EnsemblGene(ensembl_stable_id=stable_id, species_id=species_id,import_time_stamp=db.sql.functions.now())
        return ensembl_gene_value


    def build_output_hash(self):
        lines_list = []
        entry_line_dict = {
                "patho_ensembl_gene_stable_id": self.param("patho_ensembl_gene_stable_id"),
                "host_ensembl_gene_stable_id": self.param("host_ensembl_gene_stable_id"),
                }
        lines_list.append(entry_line_dict)
        return lines_list

    def write_output(self):
        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            for entry in entries_list:
                self.dataflow(entry, 1)
        else:
            print(f"{phi_id} written to FailedJob")
            self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))
    
    def check_param(self, param):
        try:
            self.param_required(param)
        except:
            error_msg = self.param('PHI_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print(error_msg)
