
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
import csv
import codecs
import eHive

class FailedEntries(eHive.BaseRunnable):
    """captures failed jobs and handles reporting"""


    def fetch_input(self):
        self.warning("FailedEntries")

    def run(self):
        phi_id = self.param("PHI_id")
        output_text = "/hps/software/users/ensembl/repositories/microbes/mcarbajo/checkout/ensembl-production-imported/lib/python/ensembl/microbes/auxiliary_files/PHIbase2/failed_entries.txt"
        msg = self.param('uncomplete_entry')
        with open(output_text, 'a+') as f:
            f.write(msg + "\n")
            f.close()
