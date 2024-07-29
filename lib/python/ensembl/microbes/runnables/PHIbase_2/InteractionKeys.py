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
import sqlalchemy as db
import sqlalchemy_utils as db_utils
import pymysql
import eHive
import ensembl.microbes.auxiliary_files.PHIbase2.interaction_DB_models as interaction_db_models
import ensembl.microbes.auxiliary_files.PHIbase2.ColumnMapper as col_map
import csv

from xml.etree import ElementTree

from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

pymysql.install_as_MySQLdb()

class InteractionKeys(eHive.BaseRunnable):
    """InteractionKeys writer to mysql DB"""

    def param_defaults(self):
        return {
            'inputfile'     : '#inputfile#',
            'registry'      : '#registry#',
            'source_db'     : '#source_db#',
            'fails_folder'  : '#fails_folder#',
        }

    def fetch_input(self):
        self.warning("Fetch InteractionKeys")
           
    
    def run(self):
        self.warning("write Keys run")
        self.param('entries_to_delete',{})
        self.insert_key_values()

    def insert_key_values(self):
        
        p2p_db_url, ncbi_tax_url, meta_db_url = self.read_registry()
        
        engine = db.create_engine(p2p_db_url)
        Session = sessionmaker(bind=engine, autoflush=False)
        session = Session()
        
        source_db = self.param('source_db')
        cm = col_map.ColumnMapper(source_db)
        key_dict = cm.keys_descriptions
        key_list = []
        for k_name in key_dict:
            self.warning("adding key " + k_name )
            key_value = self.get_key_value(session, k_name, key_dict[k_name])
            session.add(key_value)
            key_list.append(k_name)

        self.param('key_list', key_list)
        try:
            session.commit()
            self.warning("session commited")
        except pymysql.err.IntegrityError as e:
            print(e)
            session.rollback()
            self.clean_entry(engine)
        except exc.IntegrityError as e:
            print(e)
            session.rollback()
            self.clean_entry(engine)
        except Exception as e:
            print(e)
            session.rollback()
            self.clean_entry(engine)


    def get_key_value(self, session, k_name, k_desc):
        try:
            key_value = session.query(interaction_db_models.MetaKey).filter_by(name=k_name).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for {k_name}")
        except NoResultFound:
            key_value = interaction_db_models.MetaKey(name=k_name, description=k_desc)
            if 'MetaKey' in self.param('entries_to_delete'):
                added_values_list = self.param('entries_to_delete')['MetaKey']
                added_values_list.append(k_name)   
                self.add_stored_value('MetaKey',added_values_list)
            else:
                self.add_stored_value('MetaKey', [k_name])
        return key_value
        
        
    def clean_entry(self, engine):
        metadata = db.MetaData()
        connection = engine.connect()
        print("CLEAN ENTRY:")
        entries_to_delete = self.param('entries_to_delete')
        print(entries_to_delete)

        if "MetaKey" in entries_to_delete:
            key_list = entries_to_delete["MetaKey"]
            meta_key = db.Table('meta_key', metadata, autoload=True, autoload_with=engine)
            for mk_name in key_list:
                stmt = db.delete(meta_key).where(meta_key.c.name == mk_name)
                connection.execute(stmt)

    def add_stored_value(self, table,  id_value):
        new_values = self.param('entries_to_delete')
        new_values[table] = id_value
        self.param('entries_to_delete',new_values)

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
        reg_file.close()
        return int_db_url, ncbi_tax_url, meta_db_url

    def write_output(self):
        try:
            source_db = self.param('source_db')
            cm = col_map.ColumnMapper(source_db)

            obo_file = cm.ontology_file
            self.dataflow({
                "key_list": self.param('key_list'),
                "inputfile": self.param('inputfile'),
                "obo_file": obo_file,
                "fails_folder": self.param('fails_folder'),
            },1)
        except Exception as e:
            self.dataflow({
                "key_list": self.param('key_list'),
                "inputfile": self.param('inputfile'),
                "fails_folder": self.param('fails_folder'), 
            },1)

