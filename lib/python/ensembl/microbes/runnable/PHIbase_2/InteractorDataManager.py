
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
import eHive

class InteractorDataManager(eHive.BaseRunnable):
    """Makes sure all fields related to the interactors are gathered before attempting to write to mysql DB"""

    def param_defaults(self):
        return { }

    def fetch_input(self):
        self.warning("Fetch InteractorDataManager")
        self.param('branch_to_flow_on_fail', -1)
        self.param('failed_job', '')
        phi_id = self.param_required('PHI_id')
        self.check_param('patho_ensembl_gene_stable_id')
        self.check_param('patho_molecular_structure')
        self.check_param('host_ensembl_gene_stable_id')
        self.check_param('host_molecular_structure')   

    def run(self):
        self.warning("InteractorDataManager run")

    def write_output(self):
        phi_id = self.param("PHI_id")
        if self.param('failed_job') == '': 
            print(f"+++ write output for {phi_id} JOB OK")
        else:
            output_hash = [{"uncomplete_entry": self.param('failed_job')} ]
            self.dataflow(output_hash, self.param('branch_to_flow_on_fail'))
            print(f"*** write output for {phi_id} JOB KO")
            return 

    def check_param(self, param):
        try:
            self.param_required(param)
        except:
            error_msg = self.param('PHI_id') + " entry doesn't have the required field " + param + " to attempt writing to the DB"
            self.param('failed_job', error_msg)
            print(error_msg)
