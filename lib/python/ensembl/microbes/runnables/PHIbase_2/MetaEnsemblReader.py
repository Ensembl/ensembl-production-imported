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
        self.warning("/n EntryLine run")
        phi_id = self.param('PHI_id')
        #First, try to map the strain taxonomy ID
        try:
            patho_taxon_id = int(self.param('patho_species_strain'))
            patho_division, patho_db_name = self.get_meta_ensembl_info(patho_taxon_id)
            patho_species_name = self.get_species_name(patho_taxon_id)
            patho_production_name = self.get_production_name('taxonomy_id',patho_taxon_id)
        except Exception as e:
            print ("PATHO strain taxonomy ID error (will not cause the job to fail): " + str(e))
            #if strain taxonomy has failed, try the species_taxon_id
            try:
                patho_taxon_id = int(self.param('patho_species_taxon_id'))
                print ("Trying patho_species_taxon_id: "+ str(patho_taxon_id))
                patho_division, patho_db_name = self.get_meta_ensembl_info(patho_taxon_id)
                patho_species_name = self.get_species_name(patho_taxon_id)
                patho_production_name = self.get_production_name('species_taxonomy_id',patho_taxon_id)
            except Exception as e:
                print(e)
                if patho_taxon_id is None:
                    patho_taxon_id = "EMPTY"
                err_msg = "Entry: " + phi_id + " has no identifiable PATHO taxonomy (" + str(patho_taxon_id) + ")"
                self.param('failed_job',err_msg)
                return 0, 0

        self.param("patho_division",patho_division)
       	self.param("patho_dbname",patho_db_name)
       	self.param("patho_species_name",patho_species_name)
        self.param("patho_production_name",patho_production_name)


      #Same for Host
        try:
            host_taxon_id = int(self.param('host_species_strain'))
            host_division, host_db_name = self.get_meta_ensembl_info(host_taxon_id)
            host_species_name = self.get_species_name(host_taxon_id)
            host_production_name = self.get_production_name('taxonomy_id',host_taxon_id)
        except Exception as e:
            print ("HOST strain taxonomy ID error (will not cause the job to fail): " + str(e))
            #if strain taxonomy has failed, try the species_taxon_id
            try:
                host_taxon_id = int(self.param('host_species_taxon_id'))
                print ("Trying host_species_taxon_id: "+ str(host_taxon_id))
                host_division, host_db_name = self.get_meta_ensembl_info(host_taxon_id)
                host_species_name = self.get_species_name(host_taxon_id)
                host_production_name = self.get_production_name('species_taxonomy_id',host_taxon_id)
            except Exception as e:
                print(e)
                if host_taxon_id is None:
                    host_taxon_id = "EMPTY"
                err_msg = "Entry: " + phi_id + " has no identifiable HOST taxonomy (" + str(host_taxon_id) + ")" 
                self.param('failed_job',err_msg)
                return 0, 0

        self.param("host_division",host_division)
       	self.param("host_dbname",host_db_name)
        self.param("host_production_name",host_production_name)
       	self.param("host_species_name",host_species_name)

    def get_meta_ensembl_info(self, species_tax_id):
        phi_id = self.param('PHI_id')
        div_sql="SELECT DISTINCT d.short_name FROM genome g JOIN organism o USING(organism_id) JOIN division d USING(division_id) WHERE species_taxonomy_id=%d"
        self.db = pymysql.connect(host=self.param('meta_host'),user=self.param('meta_user'),db='ensembl_metadata',port=self.param('meta_port'))
        self.cur = self.db.cursor()
        try:
            self.cur.execute( div_sql % species_tax_id)
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
                self.connection_close()
        
        division = None
        if self.cur.rowcount == 1:
            division = self.cur.fetchone()[0]
        elif self.cur.rowcount == 2:
            for row in self.cur:
                if division == None and row[0] != 'EV':
                    division = row[0]
        
        if division == None:
            return 0, 0

        if division == 'EF':
            division='fungi'
        elif division == 'EV':
            division='vertebrates'
        elif division == 'EPl':
            division='plants'
        elif division == 'EPr':
            division='protists'
        elif division == 'EB':
            division='bacteria'
        elif division == 'EM':
            division='metazoa'

        core_db_name = None
        core_db_name_sql = "select gd.dbname from organism o join genome g using(organism_id) join genome_database gd using(genome_id) where o.species_taxonomy_id=%d and gd.type='core' and g.data_release_id=(select dr.data_release_id from data_release dr where is_current=1)"
        
        try:
            self.cur.execute( core_db_name_sql % species_tax_id)
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
        
        for row in self.cur:
            core_db_name = row[0]
        print (f"core_DB -- {core_db_name}")
    
        self.db.close()
        print("Division:" + str(division) + ":")

        return division, core_db_name
    
    def get_species_name(self, taxon_id):
        species_name = None
        sql="SELECT name FROM ncbi_taxa_name WHERE taxon_id=%d AND name_class='scientific name'"
        print ("*** species_name SQL: " + sql + " taxon_id=" + str(taxon_id))
        self.db = pymysql.connect(host=self.param('meta_host'),user=self.param('meta_user'),db='ncbi_taxonomy',port=self.param('meta_port'))
        self.cur = self.db.cursor()
        try:
            self.cur.execute(sql % taxon_id)
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
                self.connection_close()
        for row in self.cur:
            species_name = row[0]
        self.db.close()
        print("species_name:" + species_name)
        return species_name


    def get_production_name(self, taxon_level_field,taxon_id):
        production_name = None
        sql="SELECT name FROM organism WHERE %s=%d"
        print ("*** production_name SQL: " + sql + " taxon_id=" + str(taxon_id) + " taxon_level_field=" + str(taxon_level_field))
        self.db = pymysql.connect(host=self.param('meta_host'),user=self.param('meta_user'),db='ensembl_metadata',port=self.param('meta_port'))
        self.cur = self.db.cursor()
        try:
            self.cur.execute(sql % (taxon_level_field, taxon_id))
            self.db.commit()
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
                self.connection_close()
        for row in self.cur:
            production_name = row[0]
        self.db.close()
        print("production_name:" + production_name)
        return production_name

    def build_output_hash(self):
       lines_list = []
       entry_line_dict = {
           "patho_division": self.param("patho_division"),
           "host_division": self.param("host_division"),
           "patho_species_name": self.param("patho_species_name"),
           "host_species_name": self.param("host_species_name"),
           "patho_production_name": self.param("patho_production_name"),
           "host_production_name": self.param("host_production_name"),
           "patho_core_dbname": self.param("patho_dbname"),
           "host_core_dbname": self.param("host_dbname"),
       }
       lines_list.append(entry_line_dict)
       return lines_list

    def write_output(self):
        phi_id = self.param('PHI_id')
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
