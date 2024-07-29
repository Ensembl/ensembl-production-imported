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
from xml.etree import ElementTree
import time
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy import exc

pymysql.install_as_MySQLdb()

class EnsemblCoreReader(eHive.BaseRunnable):
    """Reads from CoreDB to gather fields to attach interactor to the ensembl_gene (mostly ensembl_gene_stable_id)"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.param('failed_job', '')
        source_db = self.param("source_db_label")
        
        entry_id = self.param_required('entry_id')

        self.check_param('interactor_A_division')
        self.check_param('interactor_B_division')
        self.check_param('interactor_A_species_taxon_id')
        self.check_param('interactor_A_dbnames_set')

    def run(self):
        if self.param('interactor_A_dbnames_set') != 0:
            self.get_ensembl_values()
        
    def get_ensembl_values(self):     
        print(self.param('entry_id'))
        interactor_A_strain_taxon_id = self.get_strain_taxon('interactor_A_species_strain')
        interactor_A_species_taxon_id = int(self.param_required('interactor_A_species_taxon_id'))
        interactor_A_taxon_ref = self.param('interactor_A_taxon_ref')
        interactor_A_molecular_id = self.param('interactor_A_molecular_id')
        interactor_A_ensembl_gene_stable_id = ''
        interactor_A_staging_url = self.get_staging_url('interactor_A_division')
        interactor_A_production_name = None
        interactor_A_scientific_name = None

        self.set_server_params('staging',interactor_A_staging_url)
        
        interactor_A_dbnames_list = self.get_names_list('interactor_A_dbnames_set')
        dbc = pymysql.connect(host=self.param('st_host'), user=self.param('st_user'), port=self.param('st_port'))
        for dbn in interactor_A_dbnames_list:
            dbn = dbn.strip()
            gg = self.getGenes_byProteinID_byDB(dbc, dbn, interactor_A_molecular_id)
            for r in gg:
                if interactor_A_scientific_name == None:
                    print("-- NONE predefined interactor_A_scientific_name ; trying db: " + str(dbn))
                    interactor_A_ensembl_gene_stable_id = r[2]
                    interactor_A_production_name = r[0]
                    interactor_A_scientific_name = r[1]
                else:
                    print("*:predefined interactor_A_scientific_name :" + str(interactor_A_scientific_name))
                        
        dbc.close()
        print("interactor_A_ensembl_gene_stable_id :" + interactor_A_ensembl_gene_stable_id 
                + ": interactor_A_production_name :" + str(interactor_A_production_name)
                + ": interactor_A_scientific_name :" + str(interactor_A_scientific_name))
        interactor_B_ensembl_gene_stable_id = ''
        interactor_B_production_name = None
        interactor_B_scientific_name = None
        if self.param('interactor_B_interactor_type') == 'synthetic':
           pass
        else:    
            interactor_B_strain_taxon_id = self.get_strain_taxon('interactor_B_species_strain')
            interactor_B_species_taxon_id = int(self.param('interactor_B_species_taxon_id'))
            interactor_B_taxon_ref = self.param('interactor_B_taxon_ref')
            interactor_B_staging_url = self.get_staging_url('interactor_B_division')
            self.set_server_params('staging', interactor_B_staging_url)
            dbc = pymysql.connect(host=self.param('st_host'), user=self.param('st_user'), port=self.param('st_port'))
            try:
                interactor_B_dbnames_list = self.get_names_list('interactor_B_dbnames_set')
                interactor_B_molecular_id = self.param('interactor_B_molecular_id')
                taxon_id = self.get_taxon_id(interactor_B_taxon_ref,interactor_B_strain_taxon_id,interactor_B_species_taxon_id)
                for dbn in interactor_B_dbnames_list: 
                    dbn = dbn.strip()
                    gg = self.getGenes_byProteinID_byDB(dbc, dbn, interactor_B_molecular_id)
                    for r in gg:
                        if  interactor_B_scientific_name == None:
                            interactor_B_ensembl_gene_stable_id = r[2]
                            interactor_B_production_name = r[0]
                            interactor_B_scientific_name = r[1]
                
                if not interactor_B_ensembl_gene_stable_id:
                    interactor_B_ensembl_gene_stable_id = self.undetermined_gene_name(self.param("interactor_B_origin_name"))

            except Exception as e:
                print(e)
                interactor_B_ensembl_gene_stable_id = self.undetermined_gene_name(self.param("interactor_B_origin_name"))
                print("interactor_B_ensembl_gene_stable_id = " + interactor_B_ensembl_gene_stable_id)
            dbc.close()        
        if interactor_A_ensembl_gene_stable_id == '':
            error_msg = self.param('entry_id') + " entry fail. Couldn't map UniProt " + interactor_A_molecular_id + " to any Ensembl gene"
            self.param('failed_job', error_msg)
            print(error_msg)
        else:
            print("** " + interactor_A_molecular_id + " mapped to " + interactor_A_ensembl_gene_stable_id + " **")
            self.param("interactor_A_ensembl_id",interactor_A_ensembl_gene_stable_id)
            self.param("interactor_A_production_name",interactor_A_production_name)
            self.param("interactor_A_scientific_name",interactor_A_scientific_name)

       # if "UNDETERMINED" not in interactor_B_ensembl_gene_stable_id: #Unfortunate double negation. Enters only if the interactor_B_stable_id is defined
        self.param("interactor_B_production_name",interactor_B_production_name)
        self.param("interactor_B_scientific_name",interactor_B_scientific_name)
        self.param("interactor_B_ensembl_id",interactor_B_ensembl_gene_stable_id)
           
    def getGenes_byProteinID_byDB(self, db_connection, dbname: str, accession_id: str):
        results = []

        gbp_sql = ("SELECT mpn.meta_value production_name, mscn.meta_value scientific_name, g.stable_id gene_stable_id " 
                     "FROM xref x "
                       "INNER JOIN object_xref ox using(xref_id) "
                       "INNER JOIN translation tr on tr.translation_id = ox.ensembl_id "
                       "INNER JOIN transcript t using(transcript_id) "
                       "INNER JOIN seq_region sr using(seq_region_id) "
                       "INNER JOIN coord_system cs using(coord_system_id) "
                       "INNER JOIN meta mpn ON mpn.species_id=cs.species_id  "
                       "INNER JOIN gene g on g.gene_id = t.gene_id "
                       "INNER JOIN meta mscn ON mscn.species_id=cs.species_id "
                     "WHERE x.dbprimary_acc = '" + accession_id + 
                       "' AND ox.ensembl_object_type = 'Translation' "
                       "AND mpn.meta_key = 'species.production_name' " 
                       "AND mscn.meta_key='species.scientific_name' "
                       "AND x.external_db_id in "
                       "(SELECT external_db_id FROM external_db WHERE db_name LIKE '%UniProt%');")

        try:
            db_connection.select_db(dbname)
            with db_connection.cursor() as cur:
                cur.execute(gbp_sql)
                results = [r for r in cur.fetchall()]
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
                self.connection_close()
        return results

    def get_species_names_by_EnsemblID(self, db_connection, dbname: str, ensembl_id: str):
        results = []   

        gbp_sql = ("SELECT mpn.meta_value production_name, mscn.meta_value scientific_name, g.stable_id gene_stable_id FROM gene g "
                          "INNER JOIN seq_region sr using(seq_region_id) "
                          "INNER JOIN coord_system cs using(coord_system_id) "
                          "INNER JOIN meta mpn using(species_id) "
                          "INNER JOIN meta mscn  using(species_id) "
                   "WHERE  mpn.meta_key = 'species.production_name' "
                          "AND mscn.meta_key='species.scientific_name' "
                          "AND g.stable_id='" + ensembl_id + "';")

        try:
            db_connection.select_db(dbname)
            with db_connection.cursor() as cur:
                cur.execute(gbp_sql)
                results = [r for r in cur.fetchall()]
        except pymysql.Error as e:
            try:
                print ("Mysql Error:- "+str(e))
            except IndexError:
                print ("Mysql Error:- "+str(e))
                self.connection_close()
        return results

    def get_ensembl_id(self, tax_id, uniprot_id, species_production_names_list):
        #returns ensembl_id and its associated species_production_name only if a block with type 'gene' is found; 0 elsewhere
        
        for sp_prod_name in species_production_names_list:
            print("... species_prod_name:" + sp_prod_name + " uniprot_id:" + uniprot_id)
            url = "https://rest.ensembl.org/xrefs/symbol/" + sp_prod_name + "/" + uniprot_id + "?external_db=UNIPROT;content-type=application/json"
            response = requests.get(url)
            unpacked_response = response.json()
            try:
                for p in unpacked_response:
                    print("getensemblid")
                    if p['type'] == 'gene':
                        print(p['id'])
                        return p['id'], sp_prod_name
            except Exception as e:
                print("Error in response: " + str(response) + "::: " + str(e))
        return '',''

    def get_ensembl_gene_value(self, session, stable_id, species_id):
        try:
            ensembl_gene_value = session.query(interaction_db_models.EnsemblGene).filter_by(ensembl_stable_id=stable_id).one()
        except MultipleResultsFound:
            print(f"ERROR: Multiple results found for species: {species_id} - gene stable_id {stable_id}")
        except NoResultFound:
            ensembl_gene_value = interaction_db_models.EnsemblGene(ensembl_stable_id=stable_id, species_id=species_id,import_time_stamp=db.sql.functions.now())
        return ensembl_gene_value

    def get_production_names_list(self, taxon_id,dbname,staging_url):
        print("taxon_id: " + taxon_id + "\dbname: " + dbname + "\staging_url: " + staging_url ) 
        if taxon_id is None: 
            return 0 
        
        production_names_list = [] 
        sql="SELECT DISTINCT meta_value FROM meta m1 WHERE m1.meta_key='species.production_name' AND m1.species_id IN (SELECT species_id FROM meta m2 WHERE m2.meta_key like '%%taxonomy_id' and meta_value=%d);" 
        try:
            self.db = pymysql.connect(host=self.param('st_host'),user=self.param('st_user'),db=dbname,port=self.param('st_port')) 
        except pymysql.Error as e:
            #Deals with the cases of model organism DB (ie.-S.cerevisiae) which have a copy of a non vertebrate DB in the vertebrates server. 
            #After trying the new connection it reverts to the previous setting.            
            vertebrates_url = self.param("vertebrate_url")
            self.set_server_params('staging', vertebrates_url)
            self.db = pymysql.connect(host=self.param('st_host'),user=self.param('st_user'),db=dbname,port=self.param('st_port'))
            self.set_server_params('staging', staging_url)

        self.cur = self.db.cursor() 
        
        try: 
            self.cur.execute(sql % (taxon_id)) 
            self.db.commit() 
        except pymysql.Error as e: 
            try: 
                print ("Mysql Error:- "+str(e)) 
            except IndexError: 
                print ("Mysql Error:- "+str(e)) 
                self.connection_close() 
        
        for row in self.cur:
            production_names_list.append(row[0])
        self.db.close() 
        return production_names_list 

    def get_names_list(self, dbname_set):
        db_names_list = []
        try:
            dbname_set_string = self.param(dbname_set)
            unbracket_dbname_set_string = dbname_set_string[1:-1]
            db_names_list = unbracket_dbname_set_string.replace("'", "").split(',')
        except:
            raise Exception(dbname_set + " is empty")
        return db_names_list

    def get_staging_url(self, division_arg):
        division = self.param(division_arg)
        staging_url = ''
        if division == 'vertebrates':
            staging_url = self.param('vertebrate_url')
        elif division == 'bacteria':
            staging_url = self.param('bacteria_url')
        else:
            staging_url = self.param('non_vertebrate_url')
        return staging_url

    def set_server_params(self,server, url):
        jdbc_pattern = 'mysql://(.*?)@(.*?):(\d*)'
        (user,host,port) = re.compile(jdbc_pattern).findall(url)[0]
        if server == 'staging':
            self.param('st_user',user)   
            self.param('st_host',host)   
            self.param('st_port',int(port))
        elif server == 'meta_ensembl':
            self.param('meta_user',user)
            self.param('meta_host',host)
            self.param('meta_port',int(port))

    def get_strain_taxon(self, species_strain):
        try:
            sp_strain = int(self.param(species_strain))
        except Exception as e:  
            sp_strain = int()
        return sp_strain

    def get_taxon_id(self, taxon_ref,strain_taxon_id,species_taxon_id):
        print("taxon_ref:" + taxon_ref + " strain_taxon_id:" + str(strain_taxon_id) + " species_taxon_id:" + str(species_taxon_id))
        if taxon_ref == 'taxonomy_id':
            return strain_taxon_id
        elif taxon_ref == 'species_taxonomy_id':
            return species_taxon_id

    def undetermined_gene_name(self, orig_name):
        return "UNDETERMINED" + "_" + self.param("entry_id") + "_" + '_'.join(orig_name.split()[:2])

    def update_molecular_id(self):
        if self.param('interactor_B_interactor_type') == 'synthetic':
            return self.param('interactor_B_molecular_id')
        try:
            print(self.param('interactor_B_molecular_id'))
            return self.param('interactor_B_molecular_id')
        except:
            return self.param('interactor_B_ensembl_id')
            '''
            try:
                print(self.param("interactor_B_production_name"))
                return self.undetermined_gene_name(self.param("interactor_B_production_name"))
            except:
                return '''''

    def build_output_hash(self):
        lines_list = []
        self.check_param('interactor_A_production_name')
        entry_line_dict = {
                "interactor_A_ensembl_id": self.param("interactor_A_ensembl_id"),
                "interactor_A_production_name": self.param("interactor_A_production_name"),
                "interactor_A_scientific_name": self.param("interactor_A_scientific_name"),
                "interactor_B_ensembl_id": self.param("interactor_B_ensembl_id"),
                "interactor_B_molecular_id": self.update_molecular_id(),
                "interactor_B_scientific_name" : self.update_scientific_name('interactor_B_scientific_name','interactor_B_origin_name'),
                "interactor_B_production_name": self.param("interactor_B_production_name"),
                }

        lines_list.append(entry_line_dict)
        return lines_list

    def write_output(self):
        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            for entry in entries_list:
                self.dataflow(entry, 1)
            print("-----------------------   CoreReader " + self.param('entry_id') + "  ----------------------------")
        else:
            print(self.param('entry_id') + " written to FailedJob")
            self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))
            print("*******************       FAILED JOB  " + self.param('entry_id') + " ***************************")
    
    def update_scientific_name(self,param, alt_param):
        try:
            if self.param(param) is None or self.param(param) == '':
                print(" scientific name via " + alt_param + " :"+  self.param(alt_param)) 
                return self.param(alt_param)
            print(" scientific name via TRY " + param + " :"+  self.param(param))
            return self.param(param)
        except:
            print(" scientific name via EXCEPT " + param + " :"+  self.param(param))
            return self.param(alt_param)
    
    def check_param(self, param):
        try:
            if self.param('interactor_A_dbnames_set') == 0:
                error_msg = self.param('entry_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
                self.param('failed_job', error_msg)
                print(error_msg)
            self.param_required(param)    
        except:
            error_msg = self.param('entry_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print(error_msg)
