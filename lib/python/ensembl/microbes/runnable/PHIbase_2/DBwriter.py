
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

from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

pymysql.install_as_MySQLdb()

class DBwriter(eHive.BaseRunnable):
    """PHI-base entry writer to mysql DB"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        phi_id = self.param_required('PHI_id')
        self.warning("Fetch WRAPPED dbWriter!" + phi_id)
        self.param('branch_to_flow_on_fail', -1)
        self.param('failed_job', '')
        
        p2p_db_url = self.param_required('interactions_db_url')
        jdbc_pattern = 'mysql://(.*?):(.*?)@(.*?):(\d*)/(.*)'
        (i_user,i_pwd,i_host,i_port,i_db) = re.compile(jdbc_pattern).findall(p2p_db_url)[0]
        self.param('p2p_user',i_user)
        self.param('p2p_pwd',i_pwd)
        self.param('p2p_host',i_host)
        self.param('p2p_port',int(i_port))
        self.param('p2p_db',i_db)
           
        self.check_param('source_db_label')
        self.check_param('patho_species_taxon_id')
        self.check_param('patho_species_name')
        self.check_param('patho_division')
        self.check_param('patho_ensembl_gene_stable_id')
        self.check_param('patho_molecula_structure')

        self.check_param('host_species_taxon_id')
        self.check_param('host_species_name')
        self.check_param('host_division')
        self.check_param('host_ensembl_gene_stable_id')
        self.check_param('host_molecular_structure')

    def run(self):
        self.param('entries_to_delete',{})
        self.warning("DBWriter run")
        self.insert_new_value()
    
    def insert_new_value(self):
        p2p_db_url = self.param_required('interactions_db_url')
        
        engine = db.create_engine(p2p_db_url)
        Session = sessionmaker(bind=engine)
        session = Session()
        session2 = Session()

        phi_id = self.param('PHI_id')
        patho_species_taxon_id = int(self.param('patho_species_taxon_id'))
        patho_species_name = self.param('patho_species_name')
        patho_division = self.param('patho_division')
        host_species_taxon_id = int(self.param('host_species_taxon_id'))
        host_species_name = self.param('host_species_name')
        host_division = self.param('host_division')

        source_db_label = self.param('source_db_label')
        print(f" PHI_id = {phi_id} :: patho_species_name {patho_species_name} :: host_species_name {host_species_name} :: source_db_label {source_db_label}")
        patho_ensembl_gene_stable_id = self.param('patho_ensembl_gene_stable_id')
        host_ensembl_gene_stable_id = self.param('host_ensembl_gene_stable_id')
        patho_molecular_structure = self.param('patho_molecular_structure')
        host_molecular_structure = self.param('host_molecular_structure')

        source_db_value = self.get_source_db_value(session, source_db_label)
        pathogen_species_value = self.get_species_value(session, patho_species_taxon_id, patho_division, patho_species_name)
        host_species_value = self.get_species_value(session, host_species_taxon_id, host_division, host_species_name)

        try:
            session.add(source_db_value)
            source_db_id = source_db_value.source_db_id
            session.add(pathogen_species_value)
            session.add(host_species_value)
            session.flush()

            pathogen_species_id = pathogen_species_value.species_id
            print(f"pathogen_species_id post flush = {pathogen_species_id}")
            host_species_id = host_species_value.species_id
            print(f"host_species_id = {host_species_id}")
            patho_ensembl_gene_value = self.get_ensembl_gene_value(session2, patho_ensembl_gene_stable_id, pathogen_species_value.species_id)
            host_ensembl_gene_value = self.get_ensembl_gene_value(session2, host_ensembl_gene_stable_id, host_species_value.species_id)
            print(f"patho_ensembl_gene_stable_id {patho_ensembl_gene_stable_id}")
            print(f"host_ensembl_gene_stable_id {host_ensembl_gene_stable_id}")
            if patho_ensembl_gene_stable_id == 'UMAG_05731':
                patho_ensembl_gene_value.gene_id = 1
                #raise ValueError('A very specific bad thing happened.')
            session2.add(patho_ensembl_gene_value)
            session2.add(host_ensembl_gene_value)
            session2.commit()
        except pymysql.err.IntegrityError as e:
            print(e)
            session.rollback()
            session2.rollback()
        except exc.IntegrityError as e:
            print(e)
            session.rollback()
            session2.rollback()
        except Exception as e:
            print(e)
            session.rollback()
            session2.rollback()
        self.clean_entry(session)

    def clean_entry(self, session):
        print("CLEAN ENTRY:")
        print(self.param('entries_to_delete'))

    def add_stored_value(self, table,  id_value):
        new_values = self.param('entries_to_delete')
        new_values[table] = id_value
        self.param('entries_to_delete',new_values)

    def get_source_db_value(self, session, db_label):
        source_db_value = None
        try:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label).one()
            print(f" db_value already exists with {db_label}")  
        except MultipleResultsFound:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label).first()
            print(f" multiple db_value exist with {db_label}") 
        except NoResultFound:
            source_db_value = interaction_db_models.SourceDb(label='PHI-base', external_db='Pathogen-Host Interactions Database that catalogues experimentally verified pathogenicity.')
            print(f" A new db_value has been created with {db_label}")
            self.add_stored_value('SourceDb', [db_label])
        return source_db_value


    def get_ensembl_gene_value(self, session, stable_id, species_id):
        try:
            ensembl_gene_value = session.query(interaction_db_models.EnsemblGene).filter_by(ensembl_stable_id=stable_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for species: {species_id} - gene stable_id {stable_id}")
        except NoResultFound:
            ensembl_gene_value = interaction_db_models.EnsemblGene(ensembl_stable_id=stable_id, species_id=species_id,import_time_stamp=db.sql.functions.now())
            if 'EnsemblGene' in self.param('entries_to_delete'):
                added_values_list = self.param('entries_to_delete')['EnsemblGene']
                added_values_list.append(stable_id)
                self.add_stored_value('EnsemblGene',added_values_list)
            else:
                self.add_stored_value('EnsemblGene', [stable_id])           
        return ensembl_gene_value

    def get_species_value(self, session, species_tax_id, division, species_name):
        
        try:
            species_value = session.query(interaction_db_models.Species).filter_by(species_taxon_id=species_tax_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for {species_name} - tx {species_tax_id}")
        except NoResultFound:
            species_value = interaction_db_models.Species(ensembl_division=division, species_production_name=species_name,species_taxon_id=species_tax_id)
            if 'Species' in self.param('entries_to_delete'):
                added_values_list = self.param('entries_to_delete')['Species']
                added_values_list.append(species_tax_id)
                self.add_stored_value('Species',added_values_list)
            else:
                self.add_stored_value('Species', [species_tax_id])   
        return species_value

    #def write_output(self):

    def check_param(self, param):
        try:
            self.param_required(param)
        except:
            error_msg = self.param('PHI_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print(error_msg)
