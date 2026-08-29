# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
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
#!/usr/bin/bash

## move to appropriate dir
cd $PWD

## Source private setup:
source ./ensembl-production-metazoa-private/conf/_assembly_status.conf

# #Source ENV to run release RR tracker
source $ASM_VENV

# Codon Cluster path
export PATH=$PATH:$MYSQL_SHORT_CMDS

## Recipient email list
ENS_METAZOA_CONTACTS="lcampbell@ebi.ac.uk ensembl-metazoa@ebi.ac.uk sdyer@ebi.ac.uk"
ENS_VERTS_CONTACTS="lcampbell@ebi.ac.uk leanne@ebi.ac.uk ensembl-genebuild@ebi.ac.uk"
ENS_PLANTS_CONTACTS="lcampbell@ebi.ac.uk gnaamati@ebi.ac.uk shradha@ebi.ac.uk sdyer@ebi.ac.uk jalvarez@ebi.ac.uk"
ENS_FUNGI_CONTACTS="lcampbell@ebi.ac.uk nishadi@ebi.ac.uk mcarbajo@ebi.ac.uk"
ENS_PROTISTS_CONTACTS="lcampbell@ebi.ac.uk nishadi@ebi.ac.uk mcarbajo@ebi.ac.uk"

## Testing / Debugging
# TEMP_EMAIL="lcampbell@ebi.ac.uk"
# ENS_METAZOA_CONTACTS=$TEMP_EMAIL
# ENS_VERTS_CONTACTS=$TEMP_EMAIL
# ENS_PLANTS_CONTACTS=$TEMP_EMAIL
# ENS_FUNGI_CONTACTS=$TEMP_EMAIL
# ENS_PROTISTS_CONTACTS=$TEMP_EMAIL

# Set processing date info
DATE_RUN=$(date "+%d-%m-%Y")
# DATE_RUN="05-07-2024"
TIME=$(date "+%T")
EXEC_DIR=`readlink -f $PWD`
TRACKING_DIR="${EXEC_DIR}/asm_tracking-${DATE_RUN}"; cd $TRACKING_DIR

for STATUS_REPORT in $(ls -1 ${TRACKING_DIR}/Ensembl*_assembly_status.tsv)
do
  
  OUT_PREFIX_DIV=`echo "$STATUS_REPORT" | awk -F "/" {'print $NF'} | sed 's/_assembly_status.tsv//g'`
  
  #Generate summarised report file:
  sh assembly_status_report_parser.sh $STATUS_REPORT $OUT_PREFIX_DIV

  #Email contents of summary and attached report files to proper recipient:
  if [[ $OUT_PREFIX_DIV == "EnsemblMetazoa" ]]; then
    mailx -s "Assembly status report [ $OUT_PREFIX_DIV / cores: staging-3 ] - $DATE_RUN" $(printf -- "-a %s " EnsemblMetazoa_*) $ENS_METAZOA_CONTACTS < ${TRACKING_DIR}/EnsemblMetazoa-SummarisedReport.txt
  elif [[ $OUT_PREFIX_DIV == "EnsemblVertebrates" ]]; then
    mailx -s "Assembly status report [ $OUT_PREFIX_DIV / cores: staging-1 ] - $DATE_RUN" $(printf -- "-a %s " EnsemblVertebrates_*) $ENS_VERTS_CONTACTS < ${TRACKING_DIR}/EnsemblVertebrates-SummarisedReport.txt
  elif [[ $OUT_PREFIX_DIV == "EnsemblPlants" ]]; then
    mailx -s "Assembly status report [ $OUT_PREFIX_DIV / cores: staging-3 ] - $DATE_RUN" $(printf -- "-a %s " EnsemblPlants_*) $ENS_PLANTS_CONTACTS < ${TRACKING_DIR}/EnsemblPlants-SummarisedReport.txt
  elif [[ $OUT_PREFIX_DIV == "EnsemblFungi" ]]; then
    mailx -s "Assembly status report [ $OUT_PREFIX_DIV / cores: staging-3 ] - $DATE_RUN" $(printf -- "-a %s " EnsemblFungi_*) $ENS_FUNGI_CONTACTS < ${TRACKING_DIR}/EnsemblFungi-SummarisedReport.txt
  elif [[ $OUT_PREFIX_DIV == "EnsemblProtists" ]]; then
    mailx -s "Assembly status report [ $OUT_PREFIX_DIV / cores: staging-3 ] - $DATE_RUN" $(printf -- "-a %s " EnsemblProtists_*) $ENS_PROTISTS_CONTACTS < ${TRACKING_DIR}/EnsemblProtists-SummarisedReport.txt
  else
    echo "Unrecognised division -> '$OUT_PREFIX_DIV' - $DATE_RUN. Check here: ${TRACKING_DIR}" | mailx -s "Issue in assembly tracking $DATE_RUN." ${USER}@ebi.ac.uk
  fi
done