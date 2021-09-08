
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
import ensembl.microbes.runnable.PHIbase_2.models as models
import ensembl.microbes.runnable.PHIbase_2.interaction_DB_models as interaction_db_models
#import dbconnection
from sqlalchemy.orm import sessionmaker


pymysql.install_as_MySQLdb()

class DBwriter(eHive.BaseRunnable):
    """PHI-base entry writer to mysql DB"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.warning("Fetch WRAPPED dbWriter!")
        phi_id = self.param_required('PHI_id')
        p2p_db_url = self.param_required('interactions_db_url')
        print(f'phi_id--{phi_id}')
        print(f'p2p_db_url--{p2p_db_url}')
        jdbc_pattern = 'mysql://(.*?):(.*?)@(.*?):(\d*)/(.*)'
        (i_user,i_pwd,i_host,i_port,i_db) = re.compile(jdbc_pattern).findall(p2p_db_url)[0]
        self.param('p2p_user',i_user)
        self.param('p2p_pwd',i_pwd)
        self.param('p2p_host',i_host)
        self.param('p2p_port',int(i_port))
        self.param('p2p_db',i_db)


        core_db_url = self.param_required('core_db_url')
        print(f'core_db_url--{core_db_url}')
        (c_user,c_pwd,c_host,c_port,c_db) = re.compile(jdbc_pattern).findall(core_db_url)[0]
        self.param('core_user',c_user)
        self.param('core_pwd',c_pwd)
        self.param('core_host',c_host)
        self.param('core_port',int(c_port))
        self.param('core_db',c_db)

    def run(self):
        self.warning("DBWriter run")
        
        p2p_db_url = self.param_required('interactions_db_url')
        species_dummy_value = interaction_db_models.Species(species_id=1, ensembl_division='Protist', species_production_name='Toxoplasma gondii',species_taxon_id=5811)
        engine = db.create_engine(p2p_db_url)
        Session = sessionmaker(bind=engine)
        session = Session()
        session.add(species_dummy_value)
        species_value = session.query(interaction_db_models.Species).filter_by(species_production_name='Toxoplasma gondii').first()
        print(f'species_value AFTER ADD -- {species_value}')
        session.delete(species_dummy_value)
        species_value = session.query(interaction_db_models.Species).filter_by(species_production_name='Toxoplasma gondii').first()
        print(f'species_value AFTER DELETE -- {species_value}')


        #core_db_url = self.param_required('core_db_url')
        #meta_dummy_value = models.Meta(meta_id=408, species_id=1, meta_key='TEST_TO_DELETE', meta_value='DELETE_ME_NOW')
        #engine = db.create_engine(core_db_url)
        #Session = sessionmaker(bind=engine)
        #session = Session()
        #session.add(meta_dummy_value)
        #meta_value = session.query(models.Meta).filter_by(meta_key='TEST_TO_DELETE').first()
        #print(f'meta_value AFTER ADD -- {meta_value}')
        #session.delete(meta_dummy_value)
        #meta_value = session.query(models.Meta).filter_by(meta_key='TEST_TO_DELETE').first()
        #print(f'meta_value AFTER DELETE -- {meta_value}')
    
        #for instance in session.query(models.Meta).order_by(models.Meta.meta_id).last():
        #    print(instance.meta_id, instance.species_id, instance.meta_key, instance.meta_value)
    
    def connection_open(self):
        self.db = pymysql.connect(host=self.param("p2p_host"),user=self.param("p2p_user"),passwd=self.param("p2p_pwd"),db=self.param("p2p_db"),port=self.param("p2p_port"))
        self.cur = self.db.cursor()

    def mysql_qry(self,sql,bool): # 1 for select and 0 for insert update delete
        self.connection_open()
        try:
            self.cur.execute(sql)
            if bool:
                return self.cur.fetchall()
            else:
                self.db.commit()
            return True
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
        self.connection_close()
             
    def connection_close(self):
        self.db.close()
        
    def myql_select(self,table):
        sql =  "SELECT * FROM "+table
        return self.mysql_qry(sql,1)

    def mysql_insert(self,table,fields,values):
        sql = "INSERT INTO " + table + " (" + fields + ") VALUES (" + values + ")";
        return self.mysql_qry(sql,0)

    def mysql_update(self,table,values,conditions):
        sql = "UPDATE " + table + " SET " + values + " WHERE " + conditions
        return self.mysql_qry(sql,0)

    def mysql_delete(self,table,conditions):
        sql = "DELETE FROM " + table + " WHERE " + conditions;
        return self.mysql_qry(sql,0)
