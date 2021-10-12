
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
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

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
        jdbc_pattern = 'mysql://(.*?):(.*?)@(.*?):(\d*)/(.*)'
        (i_user,i_pwd,i_host,i_port,i_db) = re.compile(jdbc_pattern).findall(p2p_db_url)[0]
        self.param('p2p_user',i_user)
        self.param('p2p_pwd',i_pwd)
        self.param('p2p_host',i_host)
        self.param('p2p_port',int(i_port))
        self.param('p2p_db',i_db)

        core_db_url = self.param_required('core_db_url')
        (c_user,c_pwd,c_host,c_port,c_db) = re.compile(jdbc_pattern).findall(core_db_url)[0]
        self.param('core_user',c_user)
        self.param('core_pwd',c_pwd)
        self.param('core_host',c_host)
        self.param('core_port',int(c_port))
        self.param('core_db',c_db)

        meta_db_url = self.param_required('meta_ensembl_url')
        (m_user,m_pwd,m_host,m_port,m_db) = re.compile(jdbc_pattern).findall(meta_db_url)[0]
        self.param('meta_user',m_user)
        self.param('meta_host',m_host)
        self.param('meta_port',int(m_port))

    def run(self):
        self.warning("DBWriter run")
        self.insert_new_value()
    
    def insert_new_value(self):
        p2p_db_url = self.param_required('interactions_db_url')
        
        engine = db.create_engine(p2p_db_url)
        Session = sessionmaker(bind=engine)
        session = Session()
        
        phi_id = self.param_required('PHI_id')
        patho_sequence = self.param_required('patho_sequence')
        patho_gene_name = self.param_required('patho_gene_name')
        patho_species_taxon_id = int(self.param_required('patho_species_taxon_id'))
        patho_species_name = self.param_required('patho_species_name')
        patho_core_dbname = self.param_required("patho_core_dbname")
        patho_division = self.param_required("patho_division")
        host_gene_name = self.param_required('host_gene_name')
        host_species_taxon_id = int(self.param_required('host_species_taxon_id'))
        host_species_name = self.param_required('host_species_name')
        host_core_dbname = self.param_required("host_core_dbname")
        host_division = self.param_required("host_division")

        litterature_id = self.param_required('litterature_id')
        litterature_source = self.param_required('litterature_source')
        doi = self.param_required('doi')
        interaction_phenotype = self.param_required('interaction_phenotype')
        pathogen_mutant_phenotype = self.param_required('pathogen_mutant_phenotype')
        experimental_evidence = self.param_required('experimental_evidence')
        transient_assay_exp_ev = self.param_required('transient_assay_exp_ev')

        db_label = self.get_db_label()

        source_db_value = self.get_source_db_value(session, db_label)
        self.add_if_not_exists(session, source_db_value)
        
        print(f'pathogen_species_value TO ADD -- {patho_species_name} division {patho_division} coreDB {patho_core_dbname}')
        pathogen_species_value = self.get_species_value(session, patho_species_taxon_id, patho_division, patho_species_name)
        self.add_if_not_exists(session, pathogen_species_value)
        print(f'pathogen_species_value AFTER ADD -- {pathogen_species_value} with species_id {pathogen_species_value.species_id}')
        
        print(f'host_species_value TO ADD -- {host_species_name} division {host_division} coreDB {host_core_dbname}')
        host_species_value = self.get_species_value(session, host_species_taxon_id, host_division, host_species_name)
        self.add_if_not_exists(session, host_species_value)
        print(f'host_species_value AFTER ADD -- {host_species_value}  with species_id {host_species_value.species_id}')
    
        print(f'pathogen_gene_name -- {patho_gene_name} -- division -- {pathogen_species_value.ensembl_division}')
        patho_ensembl_stable_id = self.get_ensembl_id(patho_gene_name, patho_species_taxon_id, pathogen_species_value.ensembl_division)
        host_ensembl_stable_id = self.get_ensembl_id(host_gene_name, host_species_taxon_id, host_species_value.ensembl_division)

        patho_ensembl_gene_value = self.get_ensembl_gene_value(session, patho_ensembl_stable_id, pathogen_species_value.species_id)
        host_ensembl_gene_value = self.get_ensembl_gene_value(session, host_ensembl_stable_id, host_species_value.species_id)

        self.add_if_not_exists(session, patho_ensembl_gene_value)
        self.add_if_not_exists(session, host_ensembl_gene_value)

    def get_ensembl_id(self,gene_name, tax_id ,division):
        result = None
        reported = False
        gene_names = gene_name.split(';')
        for gn in gene_names:
            if " Ensembl: " in gn:
                result = gn.replace(' Ensembl: ','')
                print ('Ensembl_id:' + result)
                reported = True
                
        if not reported:
            print ("Mmmm..., we need to find out how we name this gene name in ensembl: " + gene_name)
        return result

    def add_if_not_exists(self, session, value):
        session.add(value)
        try:
            session.commit()
        except pymysql.err.IntegrityError as e:
            print(e)
            session.rollback()
        except exc.IntegrityError as e:
            print(e)
            session.rollback()
        except Exception as e:
            print(e)
            session.rollback()


    def get_db_label(self):
        return 'PHI-base'

    def get_source_db_value(self, session, db_label):
        source_db_value = None
        try:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label).one()
        except MultipleResultsFound:
            source_db_value = session.query(interaction_db_models.SourceDb).filter_by(label=db_label).first()
        except NoResultFound:
            source_db_value = interaction_db_models.SourceDb(label='PHI-base', external_db='Pathogen-Host Interactions Database that catalogues experimentally verified pathogenicity.')
        
        return source_db_value

    def get_species_value(self, session, species_tax_id, division, species_name):
        
        try:
            species_value = session.query(interaction_db_models.Species).filter_by(species_taxon_id=species_tax_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for {species_name} - tx {species_tax_id}")
        except NoResultFound:
            species_value = interaction_db_models.Species(ensembl_division=division, species_production_name=species_name,species_taxon_id=species_tax_id)
        return species_value
    

    def get_ensembl_gene_value(self, session, stable_id, species_id):
        try:
            ensembl_gene_value = session.query(interaction_db_models.EnsemblGene).filter_by(ensembl_stable_id=stable_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for species: {species_id} - gene stable_id {stable_id}")
        except NoResultFound:
            ensembl_gene_value = interaction_db_models.EnsemblGene(ensembl_stable_id=stable_id, species_id=species_id,import_time_stamp=db.sql.functions.now())
        return ensembl_gene_value


