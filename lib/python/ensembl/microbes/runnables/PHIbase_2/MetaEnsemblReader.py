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
import re
import eHive
import requests
from requests.exceptions import HTTPError
import sqlalchemy as db
import sqlalchemy_utils as db_utils
import pymysql
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

class MetaEnsemblReader(eHive.BaseRunnable):
    """Centralises querying to ensembl-metadata and post processing of the related fields"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.param('failed_job', '')
        meta_db_url = self.param_required('meta_ensembl_url')
        jdbc_pattern = 'mysql://(.*?):(.*?)@(.*?):(\d*)/(.*)'
        (m_user,m_pwd,m_host,m_port,m_db) = re.compile(jdbc_pattern).findall(meta_db_url)[0]
        self.param('meta_user',m_user)
        self.param('meta_host',m_host)
        self.param('meta_port',int(m_port))
        
        entry_id = self.param_required('entry_id')
        interactor_B_interactor_type = self.param('interactor_B_interactor_type')

        self.check_param('interactor_A_species_taxon_id')
        self.check_param('interactor_A_name')
        self.check_param('interactor_B_name')
        if interactor_B_interactor_type == 'synthetic':
            self.param('interactor_B_species_taxon_id',0)
        else:
            self.check_param('interactor_B_species_taxon_id')
        self.check_param('doi')

    def run(self):
        if self.param('failed_job') == '':
            #self.warning("\n EntryLine run")
            entry_id = self.param('entry_id')
            db_connection = pymysql.connect(host=self.param('meta_host'),user=self.param('meta_user'),db='ensembl_metadata',port=self.param('meta_port'))
            
            species_strain = self.get_strain_taxon('interactor_A_species_strain')
            species_taxon_id = int(self.param('interactor_A_species_taxon_id'))
            interactor_A_division, interactor_A_dbnames_set, interactor_A_taxon_ref = self.get_meta_values(db_connection, species_taxon_id,species_strain)        

            self.param("interactor_A_division",interactor_A_division)
       	    self.param("interactor_A_dbnames_set",interactor_A_dbnames_set)
            self.param("interactor_A_taxon_ref",interactor_A_taxon_ref) 
            
            print("running " + entry_id)
            # Same for interactor B
            interactor_B_division = ""
            interactor_B_taxon_ref = ""
            interactor_B_dbnames_set = set()
            if self.param('interactor_B_interactor_type') == 'synthetic':
                interactor_B_division = "NA"
                interactor_B_taxon_ref = "NA"
                interactor_B_dbnames_set = set()
            else:
                species_strain = self.get_strain_taxon('interactor_B_species_strain') 
                species_taxon_id = int(self.param('interactor_B_species_taxon_id'))
                interactor_B_division, interactor_B_dbnames_set, interactor_B_taxon_ref = self.get_meta_values(db_connection, species_taxon_id, species_strain)
            
            db_connection.close()
            self.param("interactor_B_division",interactor_B_division) 
            self.param("interactor_B_dbnames_set",interactor_B_dbnames_set)
            self.param("interactor_B_taxon_ref",interactor_B_taxon_ref)

    def get_meta_values(self, db_connection, species_taxon_id, species_strain):
        division = "" 
        db_list = "" 
        db_set = set()
        entry_id = self.param('entry_id')
        #First, try to map the strain taxonomy ID
        try:
            if species_strain != 0: 
                print ("Trying strain_taxon_id: "+ str(species_strain))
                division, db_set = self.get_meta_ensembl_info(db_connection, species_strain) 
                used_taxon_ref = 'taxonomy_id'

            if division == 0 or not db_set :
                raise Exception("species_strain failed")
        #try the species taxonomy ID 
        except Exception as e:
            try:
                print ("Trying species_taxon_id: "+ str(species_taxon_id))
                division, db_set = self.get_meta_ensembl_info(db_connection, species_taxon_id)
                used_taxon_ref = 'species_taxonomy_id'
            except Exception as e:
                print(e)
                err_msg = "Entry: " + entry_id + " has no identifiable taxonomy for (" + str(species_taxon_id) + ")"
                self.param('failed_job',err_msg)
                return '',set(),''

        return division, db_set, used_taxon_ref

    def get_meta_ensembl_info(self,db_connection, tax_id):
        entry_id = self.param('entry_id')
        div_sql="SELECT DISTINCT d.short_name FROM genome g JOIN organism o USING(organism_id) JOIN division d USING(division_id) WHERE species_taxonomy_id=%d"
        self.db = db_connection
        self.cur = db_connection.cursor()
        try:
            self.cur.execute( div_sql % tax_id)
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Index Error:- "+str(e))
        
        short_div= None
        if self.cur.rowcount == 1:
            short_div = self.cur.fetchone()[0]
        elif self.cur.rowcount == 2:
            for row in self.cur:
                if short_div == None and row[0] != 'EV':
                    short_div = row[0]
        
        if short_div == None:
            return 0, 0
        division = self.get_division(short_div)
        #print("Division:" + str(division) + ":") 
        
        core_db_name = None
        core_db_name_sql = "select gd.dbname from organism o join genome g using(organism_id) join genome_database gd using(genome_id) where o.species_taxonomy_id=%d and gd.type='core' and g.data_release_id=(select MAX(dr.data_release_id) from data_release dr where is_current=1)"
        try:
            self.cur.execute( core_db_name_sql %  tax_id) 
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Index Error:- "+str(e))

        core_db_set = set()
        for row in self.cur:
            core_db_name = row[0]
            core_db_set.add(core_db_name)
            #print (f"core_DB -- {core_db_name}")
    
        return division, core_db_set
    
    def get_division(self,short_name):
        if short_name == None:
            return 0

        if short_name == 'EF':
            return 'fungi'
        elif short_name == 'EV':
            return 'vertebrates'
        elif short_name == 'EPl':
            return 'plants'
        elif short_name == 'EPr':
            return 'protists'
        elif short_name == 'EB':
            return 'bacteria'
        elif short_name == 'EM':
            return 'metazoa'

    def get_strain_taxon(self, species_strain):
        try: 
            species_strain = int(self.param(species_strain)) 
        except Exception as e:   
            species_strain = int() 
        return species_strain

    def build_output_hash(self):
       lines_list = []
       entry_line_dict = {
           "interactor_A_division": self.param("interactor_A_division"),
           "interactor_B_division": self.param("interactor_B_division"),
           "interactor_A_dbnames_set": str(self.param("interactor_A_dbnames_set")),
           "interactor_B_dbnames_set": str(self.param("interactor_B_dbnames_set")),
           "interactor_A_taxon_ref": self.param("interactor_A_taxon_ref"),
           "interactor_B_taxon_ref": self.param("interactor_B_taxon_ref"),
       }
       lines_list.append(entry_line_dict)
       return lines_list

    def write_output(self):
        self.check_param('interactor_A_species_taxon_id')
        self.check_param('interactor_B_species_taxon_id')
        self.check_param("interactor_A_division")
        self.check_param("interactor_B_division")
        self.check_param("interactor_A_dbnames_set")
        self.check_param("interactor_B_dbnames_set")
        entry_id = self.param('entry_id')

        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            for entry in entries_list:
                self.dataflow(entry, 1)
            print("--------------------------------------------------------------------------------\n") 
        else:
            print(f"{entry_id} written to FailedJob")
            self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))
            print("---------------------------****************************-------------------------\n")

    def check_param(self, param):
        try:
            self.param_required(param)
            #also checks that sets are not empty
            if "dbnames_set" in param and self.param('interactor_B_interactor_type') != 'synthetic' and self.param(param) == set():
                raise Exception('dbnames_set is empty')
        except:
            error_msg = self.param('entry_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB. "
            if param in ('interactor_A_name','interactor_B_name','interactor_A_species_taxon_id', 'interactor_B_species_taxon_id'):
                error_msg = error_msg +" Main identifier missing."
            else:                
                try:
                    if param == "interactor_A_dbnames_set":
                        error_msg = error_msg + "Could not map " + self.param('interactor_A_name') + " species_taxon(" + str(self.param('interactor_A_species_taxon_id')) + ") to Ensembl"
                    if param == "interactor_B_dbnames_set":
                        error_msg = error_msg + "Could not map " + self.param('interactor_B_name') + " species_taxon(" + str(self.param('interactor_B_species_taxon_id')) + ") to Ensembl"    
                except:
                    error_msg = error_msg + " Main identifier missing."
            self.param('failed_job', error_msg)
            print(error_msg)
