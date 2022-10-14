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
\
import os
import subprocess
import eHive

class InteractorDataManager(eHive.BaseRunnable):
    """Makes sure all fields related to the interactors are gathered before attempting to write to mysql DB"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.warning("Fetch InteractorDataManager")
        self.param('failed_job', '')
        phi_id = self.param_required('PHI_id')
        self.check_param('interactor_A_ensembl_id')
        self.check_param('interactor_A_molecular_structure')
        self.check_param('interactor_B_ensembl_id')
        self.check_param('interactor_B_molecular_structure')   

    def run(self):
        self.warning("InteractorDataManager run")
        self.get_interactor_fields()
        
    def get_interactor_fields(self):
        interactor_A_curie_type = self.param("interactor_A_curie_type")
        interactor_B_curie_type = self.param("interactor_B_curie_type")
        interactor_A_curie = interactor_A_curie_type + ":" + self.param('interactor_A_molecular_id')
        interactor_B_curie = interactor_B_curie_type + ":" + self.param('interactor_B_molecular_id')
        self.param('interactor_A_curie', interactor_A_curie)
        self.param('interactor_B_curie', interactor_B_curie)


    def get_interactor_type(self):
        #TODO Properly implement this for non PHI-base interactors
        source_db = self.param('source_db_label')
        if  source_db == 'combined_PHIbase':
            return 'protein'
        else:
            raise ValueError ("Unkonwn interactor type")

    def build_output_hash(self):
        lines_list = []
        entry_line_dict = {
                "interactor_A_interactor_type": self.param("interactor_A_interactor_type"),
                "interactor_B_interactor_type": self.param("interactor_B_interactor_type"),
                "interactor_A_curie": self.param("interactor_A_curie"),
                "interactor_B_curie": self.param("interactor_B_curie"),
                }
        lines_list.append(entry_line_dict)
        return lines_list

    def write_output(self):
        phi_id = self.param('PHI_id')
        if self.param('failed_job') == '':
            entries_list = self.build_output_hash()
            for entry in entries_list:
                self.dataflow(entry, 1)
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
