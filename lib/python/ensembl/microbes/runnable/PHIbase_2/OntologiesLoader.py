#Ensembl
#Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2021] EMBL-European Bioinformatics Institute
#
#This product includes software developed at:
#
#EMBL-European Bioinformatics Institute
#Wellcome Trust Sanger Institute

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
import csv
import ColumnMapper as col_map

from xml.etree import ElementTree
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

pymysql.install_as_MySQLdb()

class OntologiesLoader(eHive.BaseRunnable):
    """OntologiesLoader writer to mysql DB"""


    def fetch_input(self):
        self.warning("Fetch OntologiesLoader")
        self.param_required('registry')

    def run(self):
        self.warning("OntologiesLoader run")
        try:
            obo_file = self.param_required('obo_file')
            print(f"obo file {obo_file}")
            self.param('entries_to_delete',{})
            obo_dict = self.read_obo_entries()
            self.insert_obo_values(obo_dict)
        except Exception as e:
            print(f"SKIPPING LOADING ONTOLOGIES. The .obo file has not been passed as a parameter. This is not a problem as long as the necessary ontologies have previously been  loaded before.")

    def insert_obo_values(self, obo_dict):
        
        p2p_db_url, ncbi_tax_url, meta_db_url = self.read_registry()
        
        engine = db.create_engine(p2p_db_url)
        Session = sessionmaker(bind=engine)
        session = Session()
        

        source_db = self.param('source_db')
        cm = col_map.ColumnMapper(source_db)
        ontology_description = cm.ontology_description
        ontology_value = self.get_ontology_value(source_db, ontology_description, session)
        session.add(ontology_value)
        session.flush()
        onto_id = ontology_value.ontology_id

        for t_id in obo_dict:
            onto_term_value = self.get_ontology_term_value(session, onto_id, t_id, obo_dict[t_id])
            session.add(onto_term_value)
    
        try:
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

    def get_ontology_value(self, source_db, o_description, session):
        ontology_value = None
        
        try:
            ontology_value = session.query(interaction_db_models.Ontology).filter_by(name=source_db).one()
            print(f" ontology_value already exists for source_db {source_db} ")
        except MultipleResultsFound:
            ontology_value = session.query(interaction_db_models.Ontology).filter_by(name=source_db).first()
            print(f" multiple ontologies for DB  {source_db}")
        except NoResultFound:
            ontology_value = interaction_db_models.Ontology(name=source_db, description=o_description)
            print(f" A new ontology_value has been created with name {source_db} and description {o_description}")
            self.add_stored_value('Ontology', [source_db])
        return ontology_value

    def get_ontology_term_value(self, session, onto_id, t_id, t_name):
        try:
            onto_term_value = session.query(interaction_db_models.OntologyTerm).filter_by(accession=t_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for {t_id}")
        except NoResultFound:
            onto_term_value = interaction_db_models.OntologyTerm(ontology_id=onto_id, accession=t_id, description=t_name)
            
            if 'OntologyTerm' in self.param('entries_to_delete'):
                added_values_list = self.param('entries_to_delete')['OntologyTerm']
                added_values_list.append({"ontology_id":onto_id, "accession":t_id})   
                self.add_stored_value('OntologyTerm',added_values_list)
            else:
                self.add_stored_value('OntologyTerm', [{"ontology_id":onto_id, "accession":t_id}])
        return onto_term_value

    def read_obo_entries(self):
        obo_file = self.param_required('obo_file')
        of = open(obo_file, 'r')        
        Lines = of.readlines()
        obo_dict = {}
        term_id = None
        term_name = None
        for line in Lines:
            if (line.startswith("[Term]")):
                term_id = None
                term_name = None
            elif (line.startswith("id:")):
                term_id = line[4:].rstrip()
            elif (line.startswith("name:")):
                term_name = line[6:].rstrip()
            elif (not line.strip()):
                if term_id is not None:
                    obo_dict[term_id] = term_name
        of.close()
        print(f"obo dict: {obo_dict}")
        return obo_dict

    def clean_entry(self, engine):
        metadata = db.MetaData()
        connection = engine.connect()
        print("CLEAN ENTRY:")
        entries_to_delete = self.param('entries_to_delete')
        print(entries_to_delete)

        if "OntologyTerm" in entries_to_delete:
            term_list = entries_to_delete["OntologyTerm"]
            ontology_term = db.Table('ontology_term', metadata, autoload=True, autoload_with=engine)
            for obo_entry in term_list:
                stmt = db.delete(ontology_term).where(ontology_term.c.ontology_id == obo_entry['ontology_id']).where(ontology_term.c.accession == obo_entry['accession'])
                connection.execute(stmt)
                print(f"term CLEANED: {obo_entry}")

        if "Ontology" in entries_to_delete:
            ontology_name = entries_to_delete["Ontology"]
            ontology = db.Table('ontology', metadata, autoload=True, autoload_with=engine)
            stmt = db.delete(ontology).where(ontology.c.name == ontology_name)
            connection.execute(stmt)
            print(f"ontology CLEANED: {ontology_name}")

    def add_stored_value(self, table,  value_dict):
        new_values = self.param('entries_to_delete')
        new_values[table] = value_dict
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
        print(f"int_db_url {int_db_url}")
        return int_db_url, ncbi_tax_url, meta_db_url
