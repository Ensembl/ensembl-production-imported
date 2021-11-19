
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
        self.check_param('doi')
    
    def run(self):
        self.param('entries_to_delete',{})
        self.warning("DBWriter run")
        self.insert_new_value()

    def insert_new_value(self):
        p2p_db_url = self.param_required('interactions_db_url')
        
        engine = db.create_engine(p2p_db_url)
        Session = sessionmaker(bind=engine)
        session = Session()

        phi_id = self.param('PHI_id')
        patho_species_taxon_id = int(self.param('patho_species_taxon_id'))
        patho_species_name = self.param('patho_species_name')
        patho_division = self.param('patho_division')
        host_species_taxon_id = int(self.param('host_species_taxon_id'))
        host_species_name = self.param('host_species_name')
        host_division = self.param('host_division')
        source_db_label = self.param('source_db_label')
        patho_ensembl_gene_stable_id = self.param('patho_ensembl_gene_stable_id')
        host_ensembl_gene_stable_id = self.param('host_ensembl_gene_stable_id')
        patho_structure = self.param('patho_molecular_structure')
        host_structure = self.param('host_molecular_structure')
        patho_type = self.param("patho_interactor_type")
        patho_curie = self.param("patho_curie")
        patho_names = self.param("patho_other_names")
        host_type = self.param("host_interactor_type")
        host_curie = self.param("host_curie")
        host_names = self.param("host_other_names")
        doi = self.param("doi")
        print(f" PHI_id = {phi_id} :: patho_species_name {patho_species_name} :: host_species_name {host_species_name} :: source_db_label {source_db_label}")
        
        source_db_value = self.get_source_db_value(session, source_db_label)
        pathogen_species_value = self.get_species_value(session, patho_species_taxon_id, patho_division, patho_species_name)
        host_species_value = self.get_species_value(session, host_species_taxon_id, host_division, host_species_name)
        
        try:
            session.add(source_db_value)
            source_db_id = source_db_value.source_db_id
            session.add(pathogen_species_value)
            session.add(host_species_value)
            session.flush()

            patho_ensembl_gene_value = self.get_ensembl_gene_value(session, patho_ensembl_gene_stable_id, pathogen_species_value.species_id)
            host_ensembl_gene_value = self.get_ensembl_gene_value(session, host_ensembl_gene_stable_id, host_species_value.species_id)
            session.add(patho_ensembl_gene_value)
            session.add(host_ensembl_gene_value)
            session.flush()
                
            pathogen_curated_interactor = self.get_interactor_value(session, patho_type, patho_curie, patho_names, patho_structure, patho_ensembl_gene_value.gene_id)
            host_curated_interactor = self.get_interactor_value(session, host_type, host_curie, host_names, host_structure, host_ensembl_gene_value.gene_id)
            session.add(pathogen_curated_interactor)
            session.add(host_curated_interactor)
            session.flush()
            
            patho_intctr_id = pathogen_curated_interactor.curated_interactor_id
            host_intctr_id = host_curated_interactor.curated_interactor_id
            interaction_value = self.get_interaction_value(session, patho_intctr_id, host_intctr_id, doi, source_db_value.source_db_id)    
            session.add(interaction_value)

            session.commit()
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


    def get_interaction_value(self, session, patho_intctr_id, host_intctr_id, i_doi, i_source_db_id):
        interaction_value = None
        try:
            interaction_value = session.query(interaction_db_models.Interaction).filter_by(interactor_1=patho_intctr_id, interactor_2=host_intctr_id, doi=i_doi, source_db_id=i_source_db_id).one()
            print(f" interaction_value already exists with {curie}")
        except MultipleResultsFound:
            interaction_value = session.query(interaction_db_models.Interaction).filter_by(interactor_1=patho_intctr_id, interactor_2=host_intctr_id, doi=i_doi, source_db_id=i_source_db_id).first()
            print(f" multiple interaction_value exist with {curie}")
        except NoResultFound:
            interaction_value = interaction_db_models.Interaction(interactor_1=patho_intctr_id, interactor_2=host_intctr_id, doi=i_doi, source_db_id=i_source_db_id, import_timestamp=db.sql.functions.now(), meta_value_id=1)
            print(f" A new interaction_value has been created with interractors {patho_intctr_id} and {host_intctr_id} + added as stored value with {i_doi}")
            self.add_stored_value('Interaction', [{"interactor_1": patho_intctr_id, "interactor_2": host_intctr_id, "doi": i_doi, "source_db_id": i_source_db_id}])
        return interaction_value

        
    def get_interactor_value(self, session, i_type, curie, i_name, struct, gene_id):
        interactor_value = None
        try:
            interactor_value = session.query(interaction_db_models.CuratedInteractor).filter_by(curies=curie).one()
            print(f" interactor_value already exists with {curie}")
        except MultipleResultsFound:
            interactor_value = session.query(interaction_db_models.CuratedInteractor).filter_by(curies=curie).first()
        except NoResultFound:
            interactor_value = interaction_db_models.CuratedInteractor(interactor_type=i_type, curies=curie, name=i_name, molecular_structure=struct, import_timestamp=db.sql.functions.now(), ensembl_gene_id=gene_id)
            print(f" A new interactor_value has been created with {curie} + added as stored value")
            self.add_stored_value('CuratedInteractor', [curie])
        return interactor_value
        
        
    def clean_entry(self, engine):
        metadata = db.MetaData()
        connection = engine.connect()
        print("CLEAN ENTRY:")
        entries_to_delete = self.param('entries_to_delete')
        print(entries_to_delete)

        if "SourceDb" in entries_to_delete:
            source_db = db.Table('source_db', metadata, autoload=True, autoload_with=engine)
            db_label = entries_to_delete["SourceDb"]
            stmt = db.delete(source_db).where(source_db.columns.label == db_label)
            connection.execute(stmt)
            print(f"SourceDb CLEANED: {db_label}")

        if "EnsemblGene" in entries_to_delete:
            genes_list = entries_to_delete["EnsemblGene"]
            ensembl_gene = db.Table('ensembl_gene', metadata, autoload=True, autoload_with=engine)
            for stable_id in genes_list:
                stmt = db.delete(ensembl_gene).where(ensembl_gene.c.ensembl_stable_id == stable_id)
                connection.execute(stmt)
                print(f"EnsemblGene CLEANED: {stable_id}")

        if "Species" in entries_to_delete:
            species_list = entries_to_delete["Species"]
            species = db.Table('species', metadata, autoload=True, autoload_with=engine)
            for species_tax_id in species_list:
                stmt = db.delete(species).where(species.c.species_taxon_id == species_tax_id)
                connection.execute(stmt)
                print(f"Species CLEANED: {species_tax_id}")

        if "CuratedInteractor" in entries_to_delete:
            interactor_list = entries_to_delete["CuratedInteractor"]
            interactors = db.Table('curated_interactor', metadata, autoload=True, autoload_with=engine)
            for curie in interactor_list:
                stmt = db.delete(interactors).where(interactors.c.curies == curie)
                connection.execute(stmt)
                print(f"CuratedInteractor CLEANED: {curie}")
        
        if "Interaction" in entries_to_delete:
            interaction_dict = entries_to_delete["Interaction"][0]
            int_1 = interaction_dict["interactor_1"]
            int_2 = interaction_dict["interactor_2"]
            doi = interaction_dict["doi"]
            db_id = interaction_dict["source_db_id"]
            interactions = db.Table('interaction', metadata, autoload=True, autoload_with=engine)
            stmt = db.delete(interaction).where(interaction.c.interactor_1 == int_1).where(interaction.c.interactor_2 == int_2).where(interaction.c.doi == doi).where(interaction.c.source_db_id == db_id)
            connection.execute(stmt)
            print(f"Interaction CLEANED: {int_1} {int_2} {doi} {db_id}")

    def add_stored_value(self, table,  id_value):
        new_values = self.param('entries_to_delete')
        new_values[table] = id_value
        self.param('entries_to_delete',new_values)

    def get_source_db_value(self, session, db_label):
        source_db_value = None
        try:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label).one()
            print(f"db_value already exists with {db_label}")  
        except MultipleResultsFound:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label).first()
            print(f" ERROR: Multiple db_value exist with {db_label}") 
        except NoResultFound:
            source_db_value = interaction_db_models.SourceDb(label='PHI-base', external_db='Pathogen-Host Interactions Database that catalogues experimentally verified pathogenicity.')
            print(f" A new db_value has been created with {db_label} + added as stored value with sourceid {source_db_value.source_db_id}")
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
