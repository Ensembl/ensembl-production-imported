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

### Parser to summarise the results of the assembly status tracker (assembly_tracker) TSV summary report:
## Date created 28th-March-2024 (lcampbell@ebi.ac.uk)

REPORT_FILE=$1
OUTPUT_PREFIX=$2
WORK_DIR=$PWD
TIMESTAMP=`date`

if [[ -z $REPORT_FILE ]]; then
    echo "Usage: sh assembly_status_report_parser.sh <Status_report.tsv> <Output prefix name>"
    exit 1
elif [[ -z $OUTPUT_PREFIX ]]; then
    OUTPUT_PREFIX="SummarisedReport"
    echo "Output name prefix not defined. Defaulting to 'SummarisedReport'."
fi

PARSED_SUMMARY="$WORK_DIR/${OUTPUT_PREFIX}-SummarisedReport.txt"
echo -e -n "### Parsed assembly status report TSV file: '$REPORT_FILE'.\nThe following is just a summarised report; generated on $TIMESTAMP.\n\n###### Unique NCBI assembly status ######\n##columns = [assembly status, count]\n" > $PARSED_SUMMARY

ASM_QUERY_TYPE=`head -n 1 $REPORT_FILE | cut -f 1`

# Summarize non 'current' assembly status:
cut -f 10 $REPORT_FILE | grep -v -e "Asm status" | sort | uniq > status.tmp
while read STATUS
do
echo -e -n "$STATUS\t -> " >> $PARSED_SUMMARY
grep -c $STATUS $REPORT_FILE >> $PARSED_SUMMARY
if [ $STATUS != "current" ]; then 
    echo "## Following $ASM_QUERY_TYPE labelled as [$STATUS]" > ${OUTPUT_PREFIX}_Parsed_${STATUS}_cores.txt
    grep $STATUS $REPORT_FILE | cut -f1,2 >> ${OUTPUT_PREFIX}_Parsed_${STATUS}_cores.txt
fi
done < $WORK_DIR/status.tmp
rm $WORK_DIR/status.tmp


## Summarise meta data 'notes'
cut -f 11 $REPORT_FILE | grep -v -e "Asm notes" | sort | grep -v NA | uniq > asm_notes.tmp
echo -e -n "\n###### Unique assembly notes ######\n##columns = [assembly note, count]\n\n" >> $PARSED_SUMMARY
while read NOTE
do
    echo -e -n "$NOTE\t -> " >> $PARSED_SUMMARY
    grep -c -e "$NOTE" $REPORT_FILE >> $PARSED_SUMMARY
done < $WORK_DIR/asm_notes.tmp
echo -e -n "###################################\n\n" >> $PARSED_SUMMARY
rm $WORK_DIR/asm_notes.tmp

# Summarise contaminated cores:
COUNT_CONTAMINATED=`grep -e "contaminated" $REPORT_FILE | wc -l`
if [ $COUNT_CONTAMINATED -ge 1 ]; then
    echo -e -n "\n###### Following $ASM_QUERY_TYPE tagged as [contaminated] ######\n##columns = [Database | Accn, Assembly accn]\n\n" >> $PARSED_SUMMARY
    grep -e "contaminated" $REPORT_FILE | cut -f1,7 >> $PARSED_SUMMARY
    echo -e -n "###################################\n\n" >> $PARSED_SUMMARY
    echo -e -n "## Following $ASM_QUERY_TYPE tagged as [contaminated]" > ${OUTPUT_PREFIX}_Parsed_contaminated_cores.txt
    grep -e "contaminated" $REPORT_FILE | cut -f1,7 >> ${OUTPUT_PREFIX}_Parsed_contaminated_cores.txt
fi

# Understand which assemblies could be replaced entirely with different more up-to-date assemblies
COUNT_CURRENT_SUPERSEDED=`grep -e "current" $REPORT_FILE | grep -c -e "superseded by newer assembly"`
if [ $COUNT_CURRENT_SUPERSEDED -ge 1 ]; then
    echo -e -n "\n## IMPORTANT: Recovered n=$COUNT_CURRENT_SUPERSEDED 'current' assemblies but which are tagged as 'superseded' by newer assembly (diff INSDC accession):\n##columns = [Database | Accn, Assembly accn]\n\n" >> $PARSED_SUMMARY
    grep -e "current" $REPORT_FILE | grep -e "superseded by newer assembly" | cut -f1,7 >> $PARSED_SUMMARY
    echo -e -n "###################################\n\n" >> $PARSED_SUMMARY
fi

exit 0
