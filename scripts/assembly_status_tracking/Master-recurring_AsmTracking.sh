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

## Source private setup:
source ./ensembl-production-metazoa-private/conf/_assembly_status.conf

## move to appropriate dir
cd $ASM_TRACKING_WORKDIR

#Source ENV to run release RR tracker
source $ASM_VENV

# Codon Cluster path
export PATH=$PATH:$MYSQL_SHORT_CMDS

# Set processing date info
DATE_RUN=$(date "+%d-%m-%Y")
TIME=$(date "+%T")
EXEC_DIR=`readlink -f $PWD`
TRACKING_DIR="${EXEC_DIR}/asm_tracking-${DATE_RUN}"
STATUS_CHECK_HISTORY="${EXEC_DIR}/assembly_reporting_history.txt"
CURRENT_REPORT=asm_tracking-${DATE_RUN}

# Make current work dir based on date of processing
mkdir -p $TRACKING_DIR; cd $TRACKING_DIR

## Recipient email list
ENS_METAZOA_CONTACTS="lcampbell@ebi.ac.uk ensembl-metazoa@ebi.ac.uk sdyer@ebi.ac.uk"
ENS_VERTS_CONTACTS="lcampbell@ebi.ac.uk leanne@ebi.ac.uk ensembl-genebuild@ebi.ac.uk"
ENS_PLANTS_CONTACTS="lcampbell@ebi.ac.uk shradha@ebi.ac.uk sdyer@ebi.ac.uk jalvarez@ebi.ac.uk dishalodha@ebi.ac.uk"
ENS_FUNGI_CONTACTS="lcampbell@ebi.ac.uk nishadi@ebi.ac.uk mcarbajo@ebi.ac.uk"
ENS_PROTISTS_CONTACTS="lcampbell@ebi.ac.uk nishadi@ebi.ac.uk mcarbajo@ebi.ac.uk"

## Testing / Debugging
TEMP_EMAIL="lcampbell@ebi.ac.uk"
ENS_METAZOA_CONTACTS=$TEMP_EMAIL
ENS_VERTS_CONTACTS=$TEMP_EMAIL
ENS_PLANTS_CONTACTS=$TEMP_EMAIL
ENS_FUNGI_CONTACTS=$TEMP_EMAIL
ENS_PROTISTS_CONTACTS=$TEMP_EMAIL

## Gather info to decide which hosts to capture CoreDB accessions from:
# !! Staging hosts for even/odd releases sourced in private setup !!
TMP_HOST_RELEASE="$TRACKING_DIR/.release_host.tmp"
HOST_CHOICE_LOG="$EXEC_DIR/.release_host.history.txt"

# Check which ENS
echo -e -n "RELEASE_HOST\t$STAG_3\t" > $TMP_HOST_RELEASE
$STAG_3 -Ne "SHOW DATABASES LIKE 'drosophila_melanogaster%core%';" | cut -d '_' -f 5 >> $TMP_HOST_RELEASE
echo -e -n "RELEASE_HOST\t$STAG_3B\t" >> $TMP_HOST_RELEASE
$STAG_3B -Ne "SHOW DATABASES LIKE 'drosophila_melanogaster%core%';" | cut -d '_' -f 5 >> $TMP_HOST_RELEASE

ENSEMBL_VERSION_OPT_1=`head -n 1 $TMP_HOST_RELEASE | cut -f3`
ENSEMBL_VERSION_OPT_2=`tail -n 1 $TMP_HOST_RELEASE | cut -f3`

if [[ $ENSEMBL_VERSION_OPT_1 > $ENSEMBL_VERSION_OPT_2 ]]; then
	echo -e -n "Assembly tracking:$DATE_RUN $TIME\t" >> $HOST_CHOICE_LOG
	head -n 1 $TMP_HOST_RELEASE >> $HOST_CHOICE_LOG
	ENSEMBL_VERSION=$ENSEMBL_VERSION_OPT_1
	# Set specific host for both verts and nonverts
  ASMTRACK_VERT_HOST=$ASMTRACK_VERT_HOST_EVEN
	ASMTRACK_NONVERT_HOST=$ASMTRACK_NONVERT_HOST_EVEN
elif [[ $ENSEMBL_VERSION_OPT_1 < $ENSEMBL_VERSION_OPT_2 ]]; then
	echo -e -n "Assembly tracking:$DATE_RUN $TIME\t" >> $HOST_CHOICE_LOG
	tail -n 1 $TMP_HOST_RELEASE >> $HOST_CHOICE_LOG
	ENSEMBL_VERSION=$ENSEMBL_VERSION_OPT_2 ## NOT USED BUT PERHAPS USEFUL !!!
	# Set specific host for both verts and nonverts
  ASMTRACK_VERT_HOST=$ASMTRACK_VERT_HOST_ODD
	ASMTRACK_NONVERT_HOST=$ASMTRACK_NONVERT_HOST_ODD
else
  echo "Didn't detect any difference in release versions on staging hosts! Something not right. Exiting !"
  exit 0
fi


# Set specific host for both verts and nonverts - set in private conf '_assembly_status.conf'
# Now we set the connection hosts into array to pass forward
CURRENT_HOST_SERVER=($ASMTRACK_VERT_HOST $ASMTRACK_NONVERT_HOST)

# ## Gather the latest set of cores on staging host with the most current and up to date Core DBs
for HOST in ${CURRENT_HOST_SERVER[@]}
do
  CORE_LIST=all_cores_${HOST}.tmp
  CORE_LIST_DIVISION=Listed_cores_${HOST}.division.tmp

  # # ## Get full core list:
  $HOST -Ne "SHOW DATABASES LIKE '%_core_%';" | grep -v -e "collection" | grep -v -e "compara" > $CORE_LIST

  ## Get core divisions:
  while read CORE
  do
  DIVISION=`$HOST -D $CORE -Ne "SELECT meta_value FROM meta WHERE meta_key = 'species.division';"`
  echo -e -n "$CORE\t$DIVISION\n" >> $CORE_LIST_DIVISION
  done < $CORE_LIST

  ## Parse cores on division
  cut -f 2 $CORE_LIST_DIVISION | sort | uniq | grep Ensembl | xargs -I XXX sh -c "grep -e "XXX" $CORE_LIST_DIVISION | cut -f1 > XXX-${HOST}_${DATE_RUN}.corelist.txt"
done

# Delete full core and by division tmp core lists
find $TRACKING_DIR/ -type f -name "*.tmp" -delete

# Get server connection details:
ST3_PARAMS=$($ASMTRACK_NONVERT_HOST details script)
ST1_PARAMS=$($ASMTRACK_VERT_HOST details script)

#Run Assembly tracker on each set of ensembl division core list files
for CORE_LIST in Ensembl*.corelist.txt
do
DIVISION=`echo $CORE_LIST | cut -d "-" -f1`
mkdir -p $TRACKING_DIR/${DIVISION}_json_reports;
if [[ $DIVISION =~ "EnsemblVertebrates" ]]; then
  echo -e -n "Running:\n"

  #REAL CALL:
  assembly_tracker core_db --reports_dir $TRACKING_DIR/${DIVISION}_json_reports --assembly_report_name ${DIVISION}_assembly_status --input $CORE_LIST $ST1_PARAMS --log_file ${DIVISION}_assembly_tracking.log --log_file_level INFO

  ## TEST:
  # echo "assembly_tracker core_db --reports_dir $TRACKING_DIR/${DIVISION}_json_reports --assembly_report_name ${DIVISION}_assembly_status --input $CORE_LIST $ST3_PARAMS --log_file ${DIVISION}_assembly_tracking.log --log_file_level INFO"
elif [[ $DIVISION != "EnsemblVertebrates" ]]; then

  #REAL CALL:
  assembly_tracker core_db --reports_dir $TRACKING_DIR/${DIVISION}_json_reports --assembly_report_name ${DIVISION}_assembly_status --input $CORE_LIST $ST3_PARAMS --log_file ${DIVISION}_assembly_tracking.log --log_file_level INFO

  ##TEST:
  #echo "assembly_tracker core_db --reports_dir $TRACKING_DIR/${DIVISION}_json_reports --assembly_report_name ${DIVISION}_assembly_status --input $CORE_LIST $ST3_PARAMS --log_file ${DIVISION}_assembly_tracking.log --log_file_level INFO"
fi
done

## Locate set of tracking results to then summarise findings with status report parser, then tidy directory
find $TRACKING_DIR/ -type f -name "*_assembly_status.tsv" | xargs -n 1 -I XXX mv XXX ${TRACKING_DIR}/

### This is the point in which to compare between newly generated asm status report and the previous report for changes in status
CURRENT_REPORT_PRESENT=`grep -c -e "$CURRENT_REPORT" $STATUS_CHECK_HISTORY`
echo "grep -c -e $CURRENT_REPORT $STATUS_CHECK_HISTORY"
if [[ $CURRENT_REPORT_PRESENT == 0 ]]; then
    PREV_REPORT_DIR=$(tail -n 1 $STATUS_CHECK_HISTORY)
    echo $CURRENT_REPORT >> $STATUS_CHECK_HISTORY
else
    PREV_REPORT_DIR=$(grep -B 1 -e "$CURRENT_REPORT" $STATUS_CHECK_HISTORY | grep -v -e "$CURRENT_REPORT")
fi

if [[ $PREV_REPORT_DIR ]]; then
  # Loop over all assembly status TSV files
  for REPORT_TSV in $TRACKING_DIR/*_assembly_status.tsv
      do
          ASM_REPORT=$(echo $REPORT_TSV | awk -F "/" {'print $NF'})
          echo "asm report-> $ASM_REPORT"
          DIVISION=$(echo $REPORT_TSV | awk -F "/" {'print $NF'} | sed 's/_assembly_status.tsv//')
          RGREP_EXP="core_[0-9]+_[0-9].+"
          STATUS_CHANGE_REPORT="$TRACKING_DIR/${DIVISION}_status_change.txt"

          COMPARE_MESSAGE="Comparing assembly reports for status changes.\nLatest report:$CURRENT_REPORT <-> Previous report:$PREV_REPORT_DIR...\n\n"
          echo -e -n "$COMPARE_MESSAGE"
          # exit

          # Switch when no actual status change is found between reports
          NO_REAL_DIFF=1
          NO_DIFF_MESSAGE="Identical reports!\n\t >> No assembly status changes detected in latest and previous assembly status TSV. <<\n"

          echo -e -n "Division '$DIVISION':\n" > $STATUS_CHANGE_REPORT
          echo -e -n $COMPARE_MESSAGE >> $STATUS_CHANGE_REPORT

          # Check for exact match in reports, skip if identical 
          DIFFERENCE=$(diff -q $EXEC_DIR/$CURRENT_REPORT/$ASM_REPORT $EXEC_DIR/$PREV_REPORT_DIR/$ASM_REPORT)

          if [[ $DIFFERENCE ]]; then
              echo -e -n "### Latest [$CURRENT_REPORT] and previous [$PREV_REPORT_DIR] reports different.\n" >> $STATUS_CHANGE_REPORT
              echo -e -n "## Legend:\n# _Report_compared_\t_species.production_name_\t_Accession_\t_Assembly status_\t_Assembly notes_\n\n" >> $STATUS_CHANGE_REPORT

              #Loop over each entry in the CURRENT report and check for its counterpart in the previous report
              while read SINGLE_SP_STATUS
                do 
                  if [[ ! $SINGLE_SP_STATUS =~ "CoreDB" ]]; then

                    PROD_NAME=`echo "$SINGLE_SP_STATUS" | cut -d $'\t' -f1 | sed -E 's/_core.+$//g'`

                    echo -e -n "${PROD_NAME}\t" > single_cur.tmp.tsv
                    grep -E "${PROD_NAME}_${RGREP_EXP}" $EXEC_DIR/$CURRENT_REPORT/$ASM_REPORT | cut -d $'\t' -f7,10,11 >> single_cur.tmp.tsv

                    # check for same species core in current and old report to compare:
                    CORE_PRESENT=`grep -c -E "${PROD_NAME}_${RGREP_EXP}" $EXEC_DIR/$PREV_REPORT_DIR/$ASM_REPORT`

                    # If not zero we have a match and can compare report status:
                    if [[ ! $CORE_PRESENT == 0 ]]; then
                      echo -e -n "${PROD_NAME}\t" > single_past.tmp.tsv
                      grep -E "${PROD_NAME}_${RGREP_EXP}" $EXEC_DIR/$PREV_REPORT_DIR/$ASM_REPORT | cut -d $'\t' -f7,10,11 >> single_past.tmp.tsv

                      # Check if difference in status is found, then do diff if so:
                      ASM_DIFFERENCE=$(diff -q single_cur.tmp.tsv single_past.tmp.tsv)
                      if [[ $ASM_DIFFERENCE ]]; then
                        NO_REAL_DIFF=0
                        diff single_cur.tmp.tsv single_past.tmp.tsv | sed 's/</########\nLatest:/' | sed 's/>/Previous:/' | sed 's/---/ ^\n v/' >> $STATUS_CHANGE_REPORT
                        echo -e -n "\n\n" >> $STATUS_CHANGE_REPORT
                      fi
                    fi  
                  fi
                done < $EXEC_DIR/$CURRENT_REPORT/$ASM_REPORT # < This is current report as main parsing on past report

                if [[ $NO_REAL_DIFF == 1 ]]; then
                  echo -e -n "Division '$DIVISION':\n" > $STATUS_CHANGE_REPORT
                  echo -e -n "$COMPARE_MESSAGE\n$NO_DIFF_MESSAGE" >> $STATUS_CHANGE_REPORT
                fi 
                rm single_*.tmp.tsv
          else
              echo -e -n "Division '$DIVISION':\n" > $STATUS_CHANGE_REPORT
              echo -e -n "$COMPARE_MESSAGE\n$NO_DIFF_MESSAGE" >> $STATUS_CHANGE_REPORT
          fi
        echo -e -n "^^^^^^^^^^^ FINISHED COMPARISON ^^^^^^^^^^^\n\n" >> $STATUS_CHANGE_REPORT
      done
else
  # Report that there was no previous assembly tracking report located to compare the new report against
  echo -e -n "Unable to locate any assembly reports (expected to find '$PREV_REPORT_DIR') generated before today [$DATE_RUN]. Is this right ?\nNo asm status change comparison can be performed!\n\n"
fi

### Now move on to emailing reports and other files
find $TRACKING_DIR -type d -name "Ensembl*_json_reports" | xargs -n 1 -I XXX sh -c "tar -czf XXX.tar.gz XXX"
mkdir -p $TRACKING_DIR/Combined_asm_reports
find $TRACKING_DIR -type d -name "Ensembl*_json_reports" | xargs -n 1 -I XXX mv XXX $TRACKING_DIR/Combined_asm_reports/

for STATUS_REPORT in $(ls -1 ${TRACKING_DIR}/Ensembl*_assembly_status.tsv)
do

  OUT_PREFIX_DIV=`echo "$STATUS_REPORT" | awk -F "/" {'print $NF'} | sed 's/_assembly_status.tsv//g'`

  #Generate summarised report file: assembly_status_report_parser is executable
  $EXEC_DIR/assembly_status_report_parser $STATUS_REPORT $OUT_PREFIX_DIV

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
    echo "Unrecognized division -> '$OUT_PREFIX_DIV' - $DATE_RUN. Check here: ${TRACKING_DIR}" | mailx -s "Issue in assembly tracking $DATE_RUN." ${USER}@ebi.ac.uk
  fi
done

sleep 1

