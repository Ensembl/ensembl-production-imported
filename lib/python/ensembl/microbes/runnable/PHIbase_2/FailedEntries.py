
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
import csv
import codecs
import eHive

class FailedEntries(eHive.BaseRunnable):
    """captures failed jobs and handles reporting"""


    def fetch_input(self):
        self.warning("FailedEntries")

    def run(self):
        phi_id = self.param("PHI_id")
        print(f"Uncomplete entry for phi_id: {phi_id}")

