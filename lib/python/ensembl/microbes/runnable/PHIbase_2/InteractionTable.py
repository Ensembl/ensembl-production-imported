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
import eHive

class InteractionTable(eHive.BaseRunnable):
    """Manages the data that goes into the Interaction table"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.param('failed_job', '')
        phi_id = self.param_required('PHI_id')
        #self.check_param('patho_ensembl_gene_stable_id')

    def run(self):
        self.warning("InteractionTable run")
        self.get_interactor_fields()
        
    def get_interactor_fields(self):
        patho_interactor_type = self.get_interactor_type()
        host_interactor_type = self.get_interactor_type()
        patho_curie = 'prot:' + self.param('patho_uniprot_id')
        host_curie = 'prot:' + self.param('host_uniprot_id')
        self.param('patho_interactor_type', patho_interactor_type)
        self.param('host_interactor_type', host_interactor_type)
        self.param('patho_curie', patho_curie)
        self.param('host_curie', host_curie)


    def get_interactor_type(self):
        #TODO Properly implement this for non PHI-base interactors
        source_db = self.param('source_db_label')
        if  source_db == 'PHI-base':
            return 'protein'
        else:
            raise ValueError ("Unkonwn interactor type")

    def build_output_hash(self):
        lines_list = []
        entry_line_dict = {
                "patho_interactor_type": self.param("patho_interactor_type"),
                "host_interactor_type": self.param("host_interactor_type"),
                "patho_curie": self.param("patho_curie"),
                "host_curie": self.param("host_curie"),
                }
        lines_list.append(entry_line_dict)
        return lines_list

    def write_output(self):
        phi_id = self.param('PHI_id')
        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            for entry in entries_list:
                self.dataflow(entry, 1)
            print(f"{phi_id} written to DBwriter")
        else:
            #output_hash = [{"uncomplete_entry": self.param('failed_job')} ]
            print(f"{phi_id} written to FailedJob")
            self.dataflow({"uncomplete_entry": self.param('failed_job')}, self.param('branch_to_flow_on_fail'))
            return 

    def check_param(self, param):
        try:
            self.param_required(param)
        except:
            error_msg = self.param('PHI_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print(error_msg)
