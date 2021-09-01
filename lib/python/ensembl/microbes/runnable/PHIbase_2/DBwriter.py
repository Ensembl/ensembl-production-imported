
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
from datetime import datetime
import re
import ensembl.microbes.runnable.PHIbase_2.models
import dbconnection
from sqlalchemy.orm import sessionmaker


pymysql.install_as_MySQLdb()

class DBwriter(eHive.BaseRunnable):
    """PHI-base entry writer to mysql DB"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.warning("Fetch dbWriter!")
        phi_id = self.param_required('PHI_id')
        interactions_db_url = self.param_required('interactions_db_url')
        print(f'phi_id--{phi_id}')
        print(f'interactions_db_url--{interactions_db_url}')
        jdbc_pattern = 'mysql://(.*?):(.*?)@(.*?):(\d*)/(.*)'
        (user,pwd,host,port,db) = re.compile(jdbc_pattern).findall(interactions_db_url)[0]
        self.param('p2p_user',user)
        self.param('p2p_pwd',pwd)
        self.param('p2p_host',host)
        self.param('p2p_port',int(port))
        self.param('p2p_db',db)

        self.param('core_user',user)
        self.param('core_pwd',pwd)
        self.param('core_host',host)
        self.param('core_port',int(port))
        self.param('core_db',db)

    def run(self):
        self.warning("DBWriter run")
        engine = db.create_engine('mysql://ensro@mysql-ens-microbes-prod-1:4239/saccharomyces_cerevisiae_core_51_104_4')
        session = db.orm.sessionmaker(bind=engine)
        for row in session.query(Meta, Meta.meta_key).all():
            print(row.Meta, row.meta_key)

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
