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
        self.warning("EntryLine run")
	
        patho_taxon_id = int(self.param('patho_species_taxon_id'))
        try:
            patho_division, patho_db_name = self.get_meta_ensembl_info(patho_taxon_id)
            patho_species_name = self.get_species_name(patho_taxon_id)
        except Exception as e:
            print (e)
            patho_taxon_id = int(self.param('patho_species_strain'))
            patho_division, patho_db_name = self.get_meta_ensembl_info(patho_taxon_id)
            patho_species_name = self.get_species_name(patho_taxon_id)
       	
        self.param("patho_division",patho_division)
       	self.param("patho_dbname",patho_db_name)
       	self.param("patho_species_name",patho_species_name)
        
        host_taxon_id = int(self.param('host_species_taxon_id'))
        host_division, host_db_name = self.get_meta_ensembl_info(host_taxon_id)
       	host_species_name = self.get_species_name(host_taxon_id)
        
        
        self.param("host_division",host_division)
       	self.param("host_dbname",host_db_name)
       	self.param("host_species_name",host_species_name)

    def get_meta_ensembl_info(self, species_tax_id):
        div_sql="SELECT DISTINCT d.short_name FROM genome g JOIN organism o USING(organism_id) JOIN division d USING(division_id) WHERE short_name != 'EV' AND species_taxonomy_id=%d"
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
        #everything is ok, only 1 division here
            division = self.cur.fetchone()[0]
        else:
            for row in self.cur:
                if row[0] != 'EV':
                    if division == None:
                        division = row[0]
                    else:
                        # Species has more than one division
                        # remove the check for tax 559515 (uncommented for development only)
                        if species_tax_id == 559515:
                            division='EPr'
                        else:
                            print ("Division not recognised:" + str(division))
                            raise ValueError(f"Species with tax_id {species_tax_id} has more than one division")

        if division == 'EF':
            division='fungi'
        elif division == 'EPl':
            division='plants'
        elif division == 'EPr':
            division='protists'
        elif division == 'EB':
            division='bacteria'
        else:
            raise ValueError(f"Weird... That division is not supposed to be here for tax_id {species_tax_id}")

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

    def build_output_hash(self):
       lines_list = []
       entry_line_dict = {
           "patho_division": self.param("patho_division"),
           "host_division": self.param("host_division"),
           "patho_species_name": self.param("patho_species_name"),
           "host_species_name": self.param("host_species_name"),
           "patho_core_dbname": self.param("patho_dbname"),
           "host_core_dbname": self.param("host_dbname"),
       }
       lines_list.append(entry_line_dict)
       return lines_list

    def write_output(self):
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
