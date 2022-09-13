
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
        self.check_param('patho_species_taxon_id')
        self.check_param('host_species_taxon_id')
        self.check_param('patho_dbnames_set')
        self.check_param('host_dbnames_set')

    def run(self):
        self.get_values()

    def get_values(self):     
        
        patho_strain_taxon_id = self.get_strain_taxon('patho_species_strain')
        patho_species_taxon_id = int(self.param_required('patho_species_taxon_id'))
        patho_taxon_ref = self.param('patho_taxon_ref')
        patho_uniprot_id = self.param('patho_uniprot_id')

        host_strain_taxon_id = self.get_strain_taxon('host_species_strain')
        host_species_taxon_id = int(self.param_required('host_species_taxon_id'))
        host_taxon_ref = self.param('host_taxon_ref')

        patho_staging_url = self.get_staging_url('patho_division') 
        self.set_server_params('staging',patho_staging_url)
        patho_ensembl_gene_stable_id = ''
        patho_production_name = ''
        try:
            patho_ensembl_gene_stable_id = self.param_required('patho_ensembl_id')
        except Exception:
            patho_dbnames_list = self.get_names_list('patho_dbnames_set')
            taxon_id = self.get_taxon_id(patho_taxon_ref,patho_strain_taxon_id,patho_species_taxon_id)

            for dbname in patho_dbnames_list:
                if patho_ensembl_gene_stable_id == '':
                    dbname = dbname.strip()
                    print("patho db:" + dbname + ":")
                    patho_production_names_list = self.get_production_names_list(taxon_id,dbname, patho_staging_url)
                    patho_ensembl_gene_stable_id, patho_production_name = self.get_ensembl_id(taxon_id, patho_uniprot_id, patho_production_names_list)
        
        host_staging_url = self.get_staging_url('host_division')
        self.set_server_params('staging', host_staging_url)
        host_ensembl_gene_stable_id = ''
        host_production_name = ''
        try:    
            host_ensembl_gene_stable_id = self.param_required('host_ensembl_id')
        except Exception:
            print("no initial host ensembl id" )
            host_dbnames_list = self.get_names_list('host_dbnames_set')
            print("hostDB names:" + str(host_dbnames_list))
            try:
                host_uniprot_id = self.param('host_uniprot_id')
                taxon_id = self.get_taxon_id(host_taxon_ref,host_strain_taxon_id,host_species_taxon_id)
                for dbname in host_dbnames_list:
                    if not host_ensembl_gene_stable_id:
                        dbname = dbname.strip()
                        print("host db:" + dbname + ":")
                        host_production_names_list = self.get_production_names_list(taxon_id,dbname,host_staging_url)
                        host_ensembl_gene_stable_id, host_production_name = self.get_ensembl_id(taxon_id, host_uniprot_id, host_production_names_list)
                if not host_ensembl_gene_stable_id:
                    host_ensembl_gene_stable_id = "UNDETERMINED" + "_" + self.param('PHI_id') 
                    print("host_ensembl_gene_stable_id = " + host_ensembl_gene_stable_id)
            except Exception as e:
                print(e)
                host_ensembl_gene_stable_id = "UNDETERMINED" + "_" + self.param('PHI_id')
                print("host_ensembl_gene_stable_id = " + host_ensembl_gene_stable_id)

        if patho_ensembl_gene_stable_id == '':
            error_msg = self.param('PHI_id') + " entry fail. Couldn't map UniProt " + patho_uniprot_id + " to any Ensembl gene"
            self.param('failed_job', error_msg)
            print(error_msg)
        else:
            print("** " + patho_uniprot_id + " mapped to " + patho_ensembl_gene_stable_id + " **") 
            self.param("patho_ensembl_id",patho_ensembl_gene_stable_id)
            self.param("patho_species_name",patho_production_name)
        
        if "UNDETERMINED" not in host_ensembl_gene_stable_id: #Unfortunate double negation. Enters only  if the stable_id is defined
            self.param("host_species_production_name",host_production_name)
        self.param("host_ensembl_id",host_ensembl_gene_stable_id)
    

    def get_ensembl_id(self, tax_id, uniprot_id, species_production_names_list):
        #returns ensembl_id and its associated species_production_name only if a block with type 'gene' is found; 0 elsewhere
        
        for sp_prod_name in species_production_names_list:
            print("... species_prod_name:" + sp_prod_name + " uniprot_id:" + uniprot_id)
            url = "https://rest.ensembl.org/xrefs/symbol/" + sp_prod_name + "/" + uniprot_id + "?external_db=UNIPROT;content-type=application/json"
            response = requests.get(url)
            #print ("** response:" + str(response.json()))
            unpacked_response = response.json()
            try:
                for p in unpacked_response:
                    if p['type'] == 'gene':
                        print(p['id'])
                        return p['id'], sp_prod_name
            except Exception as e:
                print("Error in response: " + str(response) + "::: " + str(e))
        return '',''

    def get_ensembl_gene_value(self, session, stable_id, species_id):
        try:
            ensembl_gene_value = session.query(interaction_db_models.EnsemblGene).filter_by(ensembl_stable_id=stable_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for species: {species_id} - gene stable_id {stable_id}")
        except NoResultFound:
            ensembl_gene_value = interaction_db_models.EnsemblGene(ensembl_stable_id=stable_id, species_id=species_id,import_time_stamp=db.sql.functions.now())
        return ensembl_gene_value

    def get_production_names_list(self, taxon_id,dbname,staging_url):
        
        if taxon_id is None: 
            return 0 
        
        production_names_list = [] 
        sql="SELECT DISTINCT meta_value FROM meta m1 WHERE m1.meta_key='species.production_name' AND m1.species_id IN (SELECT species_id FROM meta m2 WHERE m2.meta_key like '%%taxonomy_id' and meta_value=%d);" 
        try:
            self.db = pymysql.connect(host=self.param('st_host'),user=self.param('st_user'),db=dbname,port=self.param('st_port')) 
        except pymysql.Error as e:
            #Deals with the cases of model organism DB (ie.-S.cerevisiae) which have a copy of a non vertebrate DB in the vertebrates server. 
            #After trying the new connection it reverts to the previous setting.            
            vertebrates_url = self.param("vertebrate_url")
            self.set_server_params('staging', vertebrates_url)
            self.db = pymysql.connect(host=self.param('st_host'),user=self.param('st_user'),db=dbname,port=self.param('st_port'))
            self.set_server_params('staging', staging_url)

        self.cur = self.db.cursor() 
        
        try: 
            self.cur.execute(sql % (taxon_id)) 
            self.db.commit() 
        except pymysql.Error as e: 
            try: 
                print ("Mysql Error:- "+str(e)) 
            except IndexError: 
                print ("Mysql Error:- "+str(e)) 
                self.connection_close() 
        
        for row in self.cur: 
            production_names_list.append(row[0])
        self.db.close() 
        return production_names_list 

    def get_names_list(self, dbname_set):
        dbname_set_string = self.param(dbname_set)
        unbracket_dbname_set_string = dbname_set_string[1:-1]
        db_names_list = unbracket_dbname_set_string.replace("'", "").split(',')
        return db_names_list

    def get_staging_url(self, division_arg):
        division = self.param(division_arg)
        staging_url = ''
        if division == 'vertebrates':
            staging_url = self.param('vertebrate_url')
        elif division == 'bacteria':
            staging_url = self.param('bacteria_url')
        else:
            staging_url = self.param('non_vertebrate_url')
        return staging_url

    def set_server_params(self,server, url):
        jdbc_pattern = 'mysql://(.*?):(.*?)@(.*?):(\d*)/(.*)'
        (user,pwd,host,port,db) = re.compile(jdbc_pattern).findall(url)[0]
        if server == 'staging':
            self.param('st_user',user)   
            self.param('st_host',host)   
            self.param('st_port',int(port))
        elif server == 'meta_ensembl':
            self.param('meta_user',user)
            self.param('meta_host',host)
            self.param('meta_port',int(port))

    def get_strain_taxon(self, species_strain):
        try:
            sp_strain = int(self.param(species_strain))
        except Exception as e:  
            sp_strain = int()
        return sp_strain

    def get_taxon_id(self, taxon_ref,strain_taxon_id,species_taxon_id):
        print("taxon_ref:" + taxon_ref + " strain_taxon_id:" + str(strain_taxon_id) + " species_taxon_id:" + str(species_taxon_id))
        if taxon_ref == 'taxonomy_id':
            return strain_taxon_id
        elif taxon_ref == 'species_taxonomy_id':
            return species_taxon_id


    def update_uniprot(self,uniprot_id):
        try:
            print ("host uniprot:" + self.param(uniprot_id))
            return self.param(uniprot_id)
        except: 
            return "UNDETERMINED"

    def update_host_species_name(self, matched_production_name, reported_name):
        try:
            print("mached_species_name:" + self.param(matched_production_name) + ":")
            return self.param(matched_production_name)
        except:
            return self.param(reported_name)

    def build_output_hash(self):
        lines_list = []
        entry_line_dict = {
                "patho_ensembl_id": self.param("patho_ensembl_id"),
                "patho_production_name": self.param("patho_species_name"),
                "host_ensembl_id": self.param("host_ensembl_id"),
                "host_uniprot_id": self.update_uniprot("host_uniprot_id"),
                "host_production_name": self.update_host_species_name("host_species_production_name","host_species_name"),
                }

        lines_list.append(entry_line_dict)
        return lines_list

    def write_output(self):
        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            for entry in entries_list:
                self.dataflow(entry, 1)
            print("-----------------------   CoreReader " + self.param('PHI_id') + "  ----------------------------")
        else:
            print(self.param('PHI_id') + " written to FailedJob")
            self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))
            print("*******************       FAILED JOB  " + self.param('PHI_id') + " ***************************")
    
    def check_param(self, param):
        try:
            self.param_required(param)
        except:
            error_msg = self.param('PHI_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print(error_msg)
