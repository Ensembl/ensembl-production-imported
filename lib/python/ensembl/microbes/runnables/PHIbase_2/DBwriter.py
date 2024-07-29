
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
import ensembl.microbes.auxiliary_files.PHIbase2.interaction_DB_models as interaction_db_models
import requests
import ensembl.microbes.auxiliary_files.PHIbase2.ColumnMapper as col_map
import csv
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
        entry_id = self.param_required('entry_id')
        self.warning("Fetch WRAPPED dbWriter!" + entry_id)
        self.param('branch_to_flow_on_fail', -1)
        self.param('failed_job', '')
        failed_entries = self.param('failed_entries',"/nfs/production/flicek/ensembl/microbes/mcarbajo/Projects/Interactions_Mapping_May_2024/failed_entries") 
        #p2p_db_url = "mysql://admin:6es4M0L1@mysql-molintxs-prod-hl.ebi.ac.uk:4716/MolIntXS_DEV"
        #self.param('interactions_db_url', p2p_db_url)

        p2p_db_url = self.param_required('interactions_db_url')
        self.param('interactions_db_url', p2p_db_url) 
        jdbc_pattern = r'mysql://(.*?):(.*?)@(.*?):(\d*)/(.*)'
        (i_user,i_pwd,i_host,i_port,i_db) = re.compile(jdbc_pattern).findall(p2p_db_url)[0]
        self.param('p2p_user',i_user)
        self.param('p2p_pwd',i_pwd)
        self.param('p2p_host',i_host)
        self.param('p2p_port',int(i_port))
        self.param('p2p_db',i_db)
           
        self.check_param('source_db_label')
        self.check_param('doi')

        self.check_param('interactor_A_species_taxon_id')
        self.check_param('interactor_A_production_name')
        self.check_param('interactor_A_division')
        self.check_param('interactor_A_ensembl_id')
        self.check_param('interactor_A_molecular_structure')

        if self.param('interactor_B_interactor_type') != 'synthetic':
            self.check_param('interactor_B_species_taxon_id')
            self.check_param('interactor_B_scientific_name')
            self.check_param('interactor_B_division')
            self.check_param('interactor_B_ensembl_id')
            self.check_param('interactor_B_molecular_structure')
        
 
    def run(self):
        self.param('entries_to_delete',{})
        self.insert_new_value()

    def insert_new_value(self):
        p2p_db_url = self.param_required('interactions_db_url')
        engine = db.create_engine(p2p_db_url)
        print ("P2P DB: " + p2p_db_url)
        Session = sessionmaker(bind=engine)
        session = Session()

        entry_id = self.param('entry_id')
        failed_entries = self.param('failed_entries')
        failed_entries_folder = str(failed_entries) + "/" + str(entry_id)
        os.makedirs(os.path.dirname(failed_entries_folder), exist_ok=True)
        print(entry_id)
        print(failed_entries_folder)
        with open(failed_entries_folder, 'w') as f:  
            with session.no_autoflush:

                f.write(entry_id)
                interactor_A_species_taxon_id = int(self.param('interactor_A_species_taxon_id'))
                interactor_A_molecular_id = self.param('interactor_A_molecular_id')
                interactor_A_production_name = self.param('interactor_A_production_name')
                interactor_A_scientific_name = self.param('interactor_A_scientific_name')
                interactor_A_division = self.param('interactor_A_division')
                interactor_A_ensembl_gene_stable_id = self.param('interactor_A_ensembl_id')
                interactor_A_structure = self.param('interactor_A_molecular_structure')
                interactor_A_interactor_type = self.param("interactor_A_interactor_type")
                interactor_A_curie = self.param("interactor_A_curie")
                interactor_A_origin_name = self.param("interactor_A_origin_name")
                interactor_B_interactor_type = self.param("interactor_B_interactor_type")
                interactor_B_molecular_id = self.param('interactor_B_molecular_id')
                interactor_B_curie = self.param("interactor_B_molecular_id")
                interactor_B_origin_name = self.param("interactor_B_origin_name")
                doi = self.param("doi")
                source_db_label = self.param('source_db_label')
                cm = col_map.ColumnMapper(source_db_label)
                source_db = self.param('source_db')
                cm = col_map.ColumnMapper(source_db)
                original_curator_db = cm.original_curator_db
                source_db_name = cm.source_db_name
                if "curated" in source_db_label:
                    source_db_label = source_db_name + "/(" + original_curator_db + " curated)"
                
                if self.param('interactor_B_interactor_type') != 'synthetic':
                    interactor_B_species_taxon_id = int(self.param('interactor_B_species_taxon_id'))
                    interactor_B_production_name = self.param('interactor_B_production_name')
                    interactor_B_scientific_name = self.param('interactor_B_scientific_name')
                    interactor_B_division = self.param('interactor_B_division')
                    interactor_B_ensembl_gene_stable_id = self.param('interactor_B_ensembl_id')
                    interactor_B_structure = self.param('interactor_B_molecular_structure')
                    interactor_B_curie = self.param("interactor_B_curie")
        
                source_db_value = self.get_source_db_value(session, source_db_label, original_curator_db)
                session.add(source_db_value)
                session.flush()
                source_db_id = source_db_value.source_db_id
                
                interactor_A_species_value = self.get_species_value(session, interactor_A_species_taxon_id, interactor_A_division, interactor_A_production_name, interactor_A_scientific_name)
                #print(entry_id + " interactor_A_production_name:" + interactor_A_production_name + " interactor_B_production_name:" + interactor_B_production_name)  
                try:
                    session.add(interactor_A_species_value)
                    if self.param('interactor_B_interactor_type') != 'synthetic':
                        interactor_B_species_value = self.get_species_value(session, interactor_B_species_taxon_id, interactor_B_division, interactor_B_production_name, interactor_B_scientific_name)
                        session.add(interactor_B_species_value)
                    session.flush()
        
                    interactor_A_ensembl_gene_value = self.get_ensembl_gene_value(session, interactor_A_ensembl_gene_stable_id, interactor_A_species_value.species_id)
                    session.add(interactor_A_ensembl_gene_value)
                    if self.param('interactor_B_interactor_type') != 'synthetic':
                        interactor_B_ensembl_gene_value = self.get_ensembl_gene_value(session, interactor_B_ensembl_gene_stable_id, interactor_B_species_value.species_id)
                        session.add(interactor_B_ensembl_gene_value)
                    session.flush()
                        
                    interactor_A_curated_interactor = self.get_interactor_value(session, interactor_A_interactor_type, interactor_A_curie, interactor_A_molecular_id, interactor_A_structure, interactor_A_ensembl_gene_value.ensembl_gene_id)
                    session.add(interactor_A_curated_interactor)
                    if self.param('interactor_B_interactor_type') != 'synthetic':
                        interactor_B_curated_interactor = self.get_interactor_value(session, interactor_B_interactor_type, interactor_B_curie, interactor_B_molecular_id, interactor_B_structure, interactor_B_ensembl_gene_value.ensembl_gene_id)
                    else:
                        interactor_B_curated_interactor = self.get_interactor_value(session, interactor_B_interactor_type, interactor_B_curie, interactor_B_origin_name, None, None)
                    session.add(interactor_B_curated_interactor)
                    session.flush()
                    
                    interactor_A_intctr_id = interactor_A_curated_interactor.curated_interactor_id
                    interactor_B_intctr_id = interactor_B_curated_interactor.curated_interactor_id
                    interaction_value = self.get_interaction_value(session, interactor_A_intctr_id, interactor_B_intctr_id, doi, source_db_value.source_db_id)    
                    session.add(interaction_value)
                    session.flush()
                    
                    interaction_id = interaction_value.interaction_id
                    key_value_pairs_dict = self.get_key_value_pairs(session)
                    
                    for key_v in key_value_pairs_dict:
                        meta_key_id = self.get_meta_key_id(session, key_v)
                        kvp_value = self.get_kv_pair_value(session, interaction_id, meta_key_id,key_value_pairs_dict[key_v])
                        session.add(kvp_value)
                    
                    session.commit()
                except pymysql.err.IntegrityError as e:
                    print(e)
                    session.rollback()
                    self.clean_entry(engine)
                    error_msg = self.param('entry_id') + " DBWriter fail: " + str(e)
                    self.param('failed_job', error_msg)
                    f.write(error_msg)
                except exc.IntegrityError as e:
                    print(e)
                    session.rollback()
                    self.clean_entry(engine)
                    error_msg = self.param('entry_id') + " DBWriter fail: " + str(e)
                    self.param('failed_job', error_msg)
                    f.write(error_msg)
                except Exception as e:
                    print(e)
                    session.rollback()
                    self.clean_entry(engine)
                    error_msg = self.param('entry_id') + " DBWriter fail: " + str(e)
                    self.param('failed_job', error_msg)
                    f.write(error_msg)

    def get_meta_key_id(self, session, key_v):
        meta_key_value = None
        try:
            meta_key_value = session.query(interaction_db_models.MetaKey).filter_by(name=key_v).one()
        except NoResultFound:
            return None
        return meta_key_value.meta_key_id	

    def get_key_value_pairs(self, session):
        key_list = self.param('key_list')
        entry_id = self.param('entry_id') 
        key_value_pairs_dict = {}
        for key in key_list:
            try:
                key_value_pairs_dict[key] = self.param_required(key)
            except Exception as e:
                pass
        return key_value_pairs_dict

    def get_kv_pair_value(self, session, int_id, key_id, mkp_value):
        kv_pair_value = None
        ontology_id = self.get_ontology_id(session)
        ontology_term_id = self.get_ontology_term_id(session, mkp_value, ontology_id)

        try:
            kv_pair_value = session.query(interaction_db_models.KeyValuePair).filter_by(interaction_id=int_id, meta_key_id=key_id, value=mkp_value).one()
        except MultipleResultsFound:
            kv_pair_value = session.query(interaction_db_models.KeyValuePair).filter_by(interaction_id=int_id, meta_key_id=key_id, value=mkp_value).first()
        except NoResultFound:
            kv_pair_value = interaction_db_models.KeyValuePair(interaction_id=int_id, meta_key_id=key_id, value=mkp_value, ontology_term_id=ontology_term_id)

            if 'KeyValuePair' in self.param('entries_to_delete'):
                added_values_list = self.param('entries_to_delete')['KeyValuePair']
                added_values_list.append({"interaction_id":int_id, "meta_key_id":key_id, "value":mkp_value, "ontology_term_id":ontology_term_id})
                self.add_stored_value('KeyValuePair',added_values_list)
            else:
                self.add_stored_value('KeyValuePair', [{"interaction_id":int_id, "meta_key_id":key_id, "value":mkp_value, "ontology_term_id":ontology_term_id}])

        return kv_pair_value
    

    def get_ontology_id(self, session):
    #TODO: Implement this for all possible ontologies
        o_name = 'PHI-DUMMY'
        ontology_value = None
        try:
            ontology_value = session.query(interaction_db_models.Ontology).filter_by(name=o_name).one()
        except NoResultFound:
            return None
        return ontology_value.ontology_id
   
    def get_ontology_term_id(self, session, o_description, o_id):
        ontology_term = None
        try:
            ontology_term = session.query(interaction_db_models.OntologyTerm).filter_by(description=o_description, ontology_id=o_id).one()
        except MultipleResultsFound:
            ontology_term = session.query(interaction_db_models.OntologyTerm).filter_by(description=o_description, ontology_id=o_id).first()
        except NoResultFound:
            return None
        
        return ontology_term.ontology_term_id

    def get_interaction_value(self, session, interactor_A_intctr_id, interactor_B_intctr_id, i_doi, i_source_db_id):
        interaction_value = None
        try:
            interaction_value = session.query(interaction_db_models.Interaction).filter_by(interactor_1=interactor_A_intctr_id, interactor_2=interactor_B_intctr_id, doi=i_doi, source_db_id=i_source_db_id).one()
        except MultipleResultsFound:
            interaction_value = session.query(interaction_db_models.Interaction).filter_by(interactor_1=interactor_A_intctr_id, interactor_2=interactor_B_intctr_id, doi=i_doi, source_db_id=i_source_db_id).first()
        except NoResultFound:
            interaction_value = interaction_db_models.Interaction(interactor_1=interactor_A_intctr_id, interactor_2=interactor_B_intctr_id, doi=i_doi, source_db_id=i_source_db_id, import_timestamp=db.sql.functions.now())
            self.add_stored_value('Interaction', [{"interactor_1": interactor_A_intctr_id, "interactor_2": interactor_B_intctr_id, "doi": i_doi, "source_db_id": i_source_db_id}])
        return interaction_value

        
    def get_interactor_value(self, session, i_type, curie, i_name, struct, gene_id):
        interactor_value = None
        try:
            interactor_value = session.query(interaction_db_models.CuratedInteractor).filter_by(curies=curie).one()
        except MultipleResultsFound:
            interactor_value = session.query(interaction_db_models.CuratedInteractor).filter_by(curies=curie).first()
        except NoResultFound:
            interactor_value = interaction_db_models.CuratedInteractor(interactor_type=i_type, curies=curie, name=i_name, molecular_structure=struct, import_timestamp=db.sql.functions.now(), ensembl_gene_id=gene_id)
            self.add_stored_value('CuratedInteractor', [curie])
        return interactor_value
        
        
    def clean_entry(self, engine):
        metadata = db.MetaData()
        connection = engine.connect()
        print("CLEAN ENTRY:")
        entries_to_delete = self.param('entries_to_delete')
        print(entries_to_delete)

        if "KeyValuePair" in entries_to_delete:
            kvp_list = entries_to_delete["KeyValuePair"]
            key_value_pair = db.Table('key_value_pair', metadata, autoload=True, autoload_with=engine)
            for kvp_dict in kvp_list:
                stmt = (db.delete(key_value_pair)
                         .where(key_value_pair.c.interaction_id == kvp_dict["interaction_id"])
                         .where(key_value_pair.c.meta_key_id == kvp_dict["meta_key_id"])
                         .where(key_value_pair.c.value == kvp_dict["value"])
                         .where(key_value_pair.c.ontology_term_id == kvp_dict["ontology_term_id"]))
                connection.execute(stmt)

        if "Interaction" in entries_to_delete:
            interaction_dict = entries_to_delete["Interaction"][0]
            int_1 = interaction_dict["interactor_1"]
            int_2 = interaction_dict["interactor_2"]
            doi = interaction_dict["doi"]
            db_id = interaction_dict["source_db_id"]
            interaction = db.Table('interaction', metadata, autoload=True, autoload_with=engine)
            stmt = db.delete(interaction).where(interaction.c.interactor_1 == int_1).where(interaction.c.interactor_2 == int_2).where(interaction.c.doi == doi).where(interaction.c.source_db_id == db_id)
            connection.execute(stmt)

        if "CuratedInteractor" in entries_to_delete:
            interactor_list = entries_to_delete["CuratedInteractor"]
            interactors = db.Table('curated_interactor', metadata, autoload=True, autoload_with=engine)
            for curie in interactor_list:
                stmt = db.delete(interactors).where(interactors.c.curies == curie)
                connection.execute(stmt)
        
            source_db = db.Table('source_db', metadata, autoload=True, autoload_with=engine)
            source_db_dict = entries_to_delete["SourceDb"][0]
            db_label = source_db_dict["db_label"]
            curator_db = source_db_dict["original_curator_db"]
            stmt = db.delete(source_db).where(source_db.columns.label == db_label).where(source_db.columns.original_curator_db == curator_db)
            connection.execute(stmt)

        if "EnsemblGene" in entries_to_delete:
            genes_list = entries_to_delete["EnsemblGene"]
            ensembl_gene = db.Table('ensembl_gene', metadata, autoload=True, autoload_with=engine)
            for stable_id in genes_list:
                stmt = db.delete(ensembl_gene).where(ensembl_gene.c.ensembl_stable_id == stable_id)
                connection.execute(stmt)

        if "Species" in entries_to_delete:
            species_list = entries_to_delete["Species"]
            species = db.Table('species', metadata, autoload=True, autoload_with=engine)
            for species_dict in species_list:
                stmt = db.delete(species).where(species.c.taxon_id == species_dict["taxon_id"]).where(species.c.scientific_name == species_dict["scientific_name"])
                connection.execute(stmt)

        if "MetaKey" in entries_to_delete:
            key_list = entries_to_delete["MetaKey"]
            meta_key = db.Table('meta_key', metadata, autoload=True, autoload_with=engine)
            for mk_name in key_list:
                stmt = db.delete(meta_key).where(meta_key.c.name == mk_name)
                connection.execute(stmt)

        if "SourceDb" in entries_to_delete:
            source_db = db.Table('source_db', metadata, autoload=True, autoload_with=engine)
            source_db_dict = entries_to_delete["SourceDb"][0]
            db_label = source_db_dict["label"]
            curator_db = source_db_dict["original_curator_db"]
            stmt = db.delete(source_db).where(source_db.columns.label == db_label).where(source_db.columns.original_curator_db == curator_db)
            connection.execute(stmt)

    def add_stored_value(self, table,  id_value):
        new_values = self.param('entries_to_delete')
        new_values[table] = id_value
        self.param('entries_to_delete',new_values)

    def get_source_db_value(self, session, db_label, curator_db):
        source_db_value = None
        try:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label, original_curator_db=curator_db).one()
        except MultipleResultsFound:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label, original_curator_db=curator_db).first()
        except NoResultFound:
            source_db_value = interaction_db_models.SourceDb(label=db_label, external_db=self.param("source_db_description"),original_curator_db=curator_db)
            self.add_stored_value('SourceDb', [{"label":db_label, "original_curator_db":curator_db}])
        return source_db_value


    def get_ensembl_gene_value(self, session, stable_id, species_id):
        try:
            ensembl_gene_value = session.query(interaction_db_models.EnsemblGene).filter_by(ensembl_stable_id=stable_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for species: {species_id} - gene stable_id {stable_id}")
        except NoResultFound:
            ensembl_gene_value = interaction_db_models.EnsemblGene(ensembl_stable_id=stable_id, species_id=species_id,import_timestamp=db.sql.functions.now())
            if 'EnsemblGene' in self.param('entries_to_delete'):
                added_values_list = self.param('entries_to_delete')['EnsemblGene']
                added_values_list.append(stable_id)
                self.add_stored_value('EnsemblGene',added_values_list)
            else:
                self.add_stored_value('EnsemblGene', [stable_id])           
        return ensembl_gene_value

    def get_species_value(self, session, species_tax_id, division, species_name, science_name):
        try:
            species_value = session.query(interaction_db_models.Species).filter_by(taxon_id=species_tax_id, scientific_name=science_name).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for {science_name} - tx {species_tax_id}")
        except NoResultFound:
            species_value = interaction_db_models.Species(ensembl_division=division, production_name=species_name, taxon_id=species_tax_id, scientific_name=science_name)
            if 'Species' in self.param('entries_to_delete'):
                added_values_list = self.param('entries_to_delete')['Species']
                added_values_list.append({"taxon_id":species_tax_id,"scientific_name":science_name})
                self.add_stored_value('Species',added_values_list)
            else:
                self.add_stored_value('Species', [{"taxon_id":species_tax_id,"scientific_name":science_name}])
        return species_value

    def write_output(self):
        if self.param('failed_job') != '':
            print(self.param('failed_job'))
            #self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))

    def check_param(self, param):
        try:
            self.param_required(param)
        except:
            error_msg = self.param('entry_id') + " DBWriter ERROR entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print("***************************************\t" + error_msg + "***************************************\t")
