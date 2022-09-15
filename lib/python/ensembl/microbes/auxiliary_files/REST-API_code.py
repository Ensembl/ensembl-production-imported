import os
import subprocess
import pymysql
import re
import requests

pymysql.install_as_MySQLdb()

#INPUT: tax_id, uniprot_id, staging_url

dbnames_list = self.get_meta_ensembl_dbnames(tax_id)

for dbname in dbnames_list:
    production_names_list = self.get_production_names_list(taxon_id,dbname, staging_url)
    ensembl_gene_stable_id, patho_production_name = self.get_ensembl_id(taxon_id, patho_uniprot_id, production_names_list)


def get_ensembl_id(self, tax_id, uniprot_id, species_production_names_list):
    for sp_prod_name in species_production_names_list:
        url = "https://rest.ensembl.org/xrefs/symbol/" + sp_prod_name + "/" + uniprot_id + "?external_db=UNIPROT;content-type=application/json"
        response = requests.get(url)
        unpacked_response = response.json()
        try:
            for p in unpacked_response:
                if p['type'] == 'gene':
                    return p['id'], sp_prod_name            
        except Exception as e:
            print(e)
    return '',''
        
        
def get_meta_ensembl_dbnames(self,tax_id):        
    core_db_name = None
    core_db_name_sql = "SELECT gd.dbname FROM organism o JOIN genome g USING(organism_id) JOIN genome_database gd USING(genome_id) WHERE o.species_taxonomy_id=%d AND gd.type='core' AND g.data_release_id=(SELECT MAX(dr.data_release_id) FROM data_release dr WHERE is_current=1)"
        try:
            self.cur.execute(core_db_name_sql %  tax_id) 
            self.db.commit()
        except pymysql.Error as e:
            print (e)
        core_db_set = set() #sets have the advantage of being unique
        for row in self.cur:
            core_db_set.add(row[0])

        self.db.close()
        return list(core_db_set) #returns set as a list
    
    
def get_production_names_list(self, taxon_id,dbname,staging_url):
    production_names_list = [] 
    sql="SELECT DISTINCT meta_value FROM meta m1 WHERE m1.meta_key='species.production_name' AND m1.species_id IN (SELECT species_id FROM meta m2 WHERE m2.meta_key like '%%taxonomy_id' and meta_value=%d);" 
    self.db = pymysql.connect(host=self.param('st_host'),user=self.param('st_user'),db=dbname,port=self.param('st_port')) 
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
