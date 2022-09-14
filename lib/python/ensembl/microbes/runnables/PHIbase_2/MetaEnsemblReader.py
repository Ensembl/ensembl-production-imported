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

        phi_id = self.param_required('PHI_id')
        self.check_param('patho_species_taxon_id')
        self.check_param('host_species_taxon_id')
        self.check_param('doi')

    def run(self):
        if self.param('failed_job') == '':
            #self.warning("\n EntryLine run")
            phi_id = self.param('PHI_id')
        
            species_strain = self.get_strain_taxon('patho_species_strain')
            species_taxon_id = int(self.param('patho_species_taxon_id'))
            patho_division, patho_db_names_set, patho_taxon_ref = self.get_meta_values(species_taxon_id,species_strain)        

            self.param("patho_division",patho_division)
       	    self.param("patho_dbnames_set",patho_db_names_set)
            self.param("patho_taxon_ref",patho_taxon_ref) 
        
            # Same for Host 
            species_strain = self.get_strain_taxon('host_species_strain') 
            species_taxon_id = int(self.param('host_species_taxon_id'))
            host_division, host_db_names_set, host_taxon_ref = self.get_meta_values(species_taxon_id, species_strain)

            self.param("host_division",host_division) 
            self.param("host_dbnames_set",host_db_names_set)
            self.param("host_taxon_ref",host_taxon_ref)

    def get_meta_values(self, species_taxon_id, species_strain):
        division = "" 
        db_list = "" 
        db_set = set()
        phi_id = self.param('PHI_id')
        #First, try to map the strain taxonomy ID
        try:
            if species_strain != 0: 
                #print ("Trying strain_taxon_id: "+ str(species_strain))
                division, db_set = self.get_meta_ensembl_info(species_strain) 
                used_taxon_ref = 'taxonomy_id'

            if division == 0 or not db_set :
                raise Exception("species_strain failed")
        #try the species taxonomy ID 
        except Exception as e:
            try:
                #print (str(e) + " Trying species_taxon_id: "+ str(species_taxon_id))
                division, db_set = self.get_meta_ensembl_info(species_taxon_id)
                used_taxon_ref = 'species_taxonomy_id'
            except Exception as e:
                print(e)
                err_msg = "Entry: " + phi_id + " has no identifiable taxonomy for (" + str(species_taxon_id) + ")"
                self.param('failed_job',err_msg)
                return 0, 0, 0

        return division, db_set, used_taxon_ref

    def get_meta_ensembl_info(self, tax_id):
        phi_id = self.param('PHI_id')
        div_sql="SELECT DISTINCT d.short_name FROM genome g JOIN organism o USING(organism_id) JOIN division d USING(division_id) WHERE species_taxonomy_id=%d"
        self.db = pymysql.connect(host=self.param('meta_host'),user=self.param('meta_user'),db='ensembl_metadata',port=self.param('meta_port'))
        self.cur = self.db.cursor()
        try:
            self.cur.execute( div_sql % tax_id)
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
                self.connection_close()
        
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
        #core_db_name_sql = "select gd.dbname from organism o  join genome g using(organism_id) join genome_database gd using(genome_id) join division d on g.division_id=d.division_id where d.short_name=%s and o.species_taxonomy_id=%d and gd.type='core' and g.data_release_id=(select MAX(dr.data_release_id) from data_release dr where is_current=1)"
        try:
            self.cur.execute( core_db_name_sql %  tax_id) 
            #self.cur.execute( core_db_name_sql % (short_div,tax_id))
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
        core_db_set = set()
        for row in self.cur:
            core_db_name = row[0]
            core_db_set.add(core_db_name)
            #print (f"core_DB -- {core_db_name}")
    
        self.db.close()
        #print("coredbset:")
        #print (core_db_set)  
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
           "patho_division": self.param("patho_division"),
           "host_division": self.param("host_division"),
           "patho_dbnames_set": str(self.param("patho_dbnames_set")),
           "host_dbnames_set": str(self.param("host_dbnames_set")),
           "patho_taxon_ref": self.param("patho_taxon_ref"),
           "host_taxon_ref": self.param("host_taxon_ref"),
       }
       lines_list.append(entry_line_dict)
       return lines_list

    def write_output(self):
        self.check_param("patho_division")
        self.check_param("host_division")
        self.check_param("patho_dbnames_set")
        self.check_param("host_dbnames_set")
        phi_id = self.param('PHI_id')

        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            for entry in entries_list:
                self.dataflow(entry, 1)
            print("--------------------------------------------------------------------------------\n") 
        else:
            print(f"{phi_id} written to FailedJob")
            self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))
            print("---------------------------****************************-------------------------\n")

    def check_param(self, param):
        try:
            self.param_required(param)
            test_param = self.param(param) 
            #also checks that sets are not empty
            if "dbnames_set" in param and not test_param:
                raise Exception('dbnames_set is empty')
        except:
            error_msg = self.param('PHI_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB. "
            if param == "patho_dbnames_set":
                error_msg = error_msg + "Could not map " + self.param('patho_species_name') + " species_taxon(" + str(self.param('patho_species_taxon_id')) + ") to Ensembl"
            if param == "host_dbnames_set":
                error_msg = error_msg + "Could not map " + self.param('host_species_name') + " species_taxon(" + str(self.param('host_species_taxon_id')) + ") to Ensembl"    
            self.param('failed_job', error_msg)
            print(error_msg)
