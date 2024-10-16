#!/usr/bin/env bash
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

#### Script to facilitate batch submission of db copy jobs on cores (host) facilitating
# core renames (target) via the production-tools 'dbcopy-client'.
# Alternativelt it facilitates listing of user specific completed copy jobs.

EMAIL_LIST="${USER}@ebi.ac.uk"  # Can be a comma separated list of email addresses the system will notify

SOURCE_HOST=$1
SOURCE=$2
TARGET_HOST=$3
TARGET=$4
OPERATION=$5
CWD=`readlink -f $PWD`
INPUT_CORES="$CWD/combined_copy_cores.tmp"
LOG_FILE="$CWD/BATCH_DB_COPY_${OPERATION^}.log"

if [[ -z $SOURCE ]] || [[ -z $SOURCE_HOST ]] || [[ -z $TARGET ]] || [[ -z $TARGET_HOST ]] || [[ -z $OPERATION ]]; then

	echo "Usage: Bulk_DB_Copy <SOURCE host> <Source Core(s) list file> <TARGET host> <Target Core(s) list file> <Operation: submit | list>"
	echo "Source (A) and Target (B) core list files should be ordered lists 'A1' => 'B1', 'A2' => 'B2' etc."
	echo "E.g. Bulk_DB_Copy me1 Source_cores.txt me2 Target_cores.txt submit"
	exit 0
fi

##Vars for STDOUT colour formatting
read -p "Do you want colourised STDOUT printing enabled ? <yes|no> " COLOUR_OPT
if [[ "$COLOUR_OPT" == "YES" ]] || [[ "$COLOUR_OPT" == "yes" ]]; then
	GREEN='\033[0;32m'
	RED='\033[0;31m'
	ORANGE='\033[0;33m'
	NC='\033[0m'
else
	GREEN=''
	RED=''
	ORANGE=''
	NC=''
fi

# Check input files passed have same length before combining
SOURCE_CORES=`readlink -f $SOURCE`
COUNT_SOURCE=`wc -l $SOURCE_CORES | cut -d " " -f1`
TARGET_CORES=`readlink -f $TARGET`
COUNT_TARGET=`wc -l $TARGET_CORES | cut -d " " -f1`

if [[ $COUNT_SOURCE -ne $COUNT_TARGET ]]; then
	echo -e -n "${RED}Number of CoreDBs in source and target core lists do not match! Ensure both lists have the same number of cores per file.${NC}\n"
	exit 1
fi

`paste $SOURCE_CORES $TARGET_CORES > $INPUT_CORES`

# Find out if we want to do a test run:
read -p "Set run as a TEST dry run ? <yes|no> " DRY_RUN
if [[ "${DRY_RUN^^}" == "YES" ]]; then
	echo -e -n "${ORANGE}SAFE RUN: DB Copy dry test run.${NC}\n\n"
elif [[ "${DRY_RUN^^}" == "NO" ]]; then
	echo -e -n "${ORANGE}REAL RUN: DB Copy will attempt to submit real jobs!!${NC}\n\n"
	sleep 2
else
	echo -e -n "${RED}Options are 'yes' or 'no'${NC}\n"
	exit 1
fi


### Main production endpoints
if [[ -e $ENSEMBL_ROOT_DIR/"ensembl-production-imported-private/lib/bash/HANDOVER.env" ]]; then
	source $ENSEMBL_ROOT_DIR/ensembl-production-imported-private/lib/bash/HANDOVER.env
	URL_SOURCE="$ENSEMBL_ROOT_DIR/ensembl-production-imported-private/lib/bash/HANDOVER.env"
elif [[ -e $CWD/"ensembl-production-imported-private/lib/bash/HANDOVER.env" ]]; then
	source $CWD/ensembl-production-imported-private/lib/bash/HANDOVER.env
	URL_SOURCE="$CWD/ensembl-production-imported-private/lib/bash/HANDOVER.env"
else
	echo -e -n "${ORANGE}Could not detect local clone of 'ensembl-production-imported-private'\n"
	echo -e -n "Attempting fresh clone now...${NC}\n\n"
	git clone -b main --depth 1 git@github.com:Ensembl/ensembl-production-imported-private.git $CWD/ensembl-production-imported-private
	source $CWD/ensembl-production-imported-private/lib/bash/HANDOVER.env
	URL_SOURCE="$CWD/ensembl-production-imported-private/lib/bash/HANDOVER.env"
fi

# Check for endpoint URL correctly defined and set for db-copy purposes
if [ $HANDOVER_BASE_URL ]; then
	echo -e -n "${GREEN}Detected dbcopy-client endpoint URL variable -> $URL_SOURCE${NC}\n"
	ENDPOINT="${HANDOVER_BASE_URL}dbcopy/requestjob"
else
	echo -e -n "${RED}Unable to detect Endpoint URL variable for dbcopy-client.\n Required repo: 'ensembl-production-imported-private'${NC}\n\n"
fi

# Check for the right modenv is located and activated
which dbcopy-client
if [ $? -eq 0 ]; then
	echo -e -n "${GREEN}ensembl_mod_env module environment detected 'ensembl/production-tools' dbcopy-client found! Proceeding...${NC}\n\n"
else
	echo -e -n "${RED}Unable to locate 'dbcopy-client' (production-tools). Ensure you have created/loaded via 'module load ensembl/production-tools'.\n"
	echo -e -n "modenv_create ensembl/production-tools production-tools && module load ensembl/production-tools'${NC}\n"
	echo -e -n "${ORANGE}See --> https://github.com/Ensembl/ensembl-mod-env/wiki/\n\n"
	exit 1
fi

# Database host expansion:
SOURCE_SERVER=`$SOURCE_HOST details host-port`
TARGET_SERVER=`$TARGET_HOST details host-port`

BREAK=0
if [[ $BREAK == 0 ]]; then
	date 2>&1 | tee $LOG_FILE
	echo -e -n "${GREEN}##Source server: $SOURCE_SERVER ----> Target server: $TARGET_SERVER${NC}\n" | tee -a $LOG_FILE
	BREAK=1
fi

while read LINE; do

	CORE=`echo "$LINE" | cut -d $'\t' -f1`
	TARGET_CORE=`echo "$LINE" | cut -d $'\t' -f2`
	
	# Only do dry run
	if [[ "${DRY_RUN^^}" == "YES" ]] && [[ "${OPERATION^^}" == "SUBMIT" ]]; then
		echo -e -n "${ORANGE}\n#Dry Run DB copy...\n\t" | tee -a $LOG_FILE
		echo -e -n "=> dbcopy-client -a $OPERATION -u $ENDPOINT --src_host $SOURCE_SERVER --tgt_host $TARGET_SERVER -e $EMAIL_LIST -r $USER -i $CORE -n $TARGET_CORE --skip-check${NC}\n\n" | tee -a $LOG_FILE
	#Make transfer with dbcopy-client
	elif [[ "${DRY_RUN^^}" == "NO" ]] && [[ "${OPERATION^^}" == "SUBMIT" ]]; then
		echo -e -n "${GREEN}#Real copy submission...\n\t"
		echo -e -n "=> dbcopy-client -a $OPERATION -u $ENDPOINT --src_host $SOURCE_SERVER --tgt_host $TARGET_SERVER -e $EMAIL_LIST -r $USER -i $CORE -n $TARGET_CORE --skip-check${NC}\n" | tee -a $LOG_FILE
		dbcopy-client -a $OPERATION -u $ENDPOINT --src_host $SOURCE_SERVER --tgt_host $TARGET_SERVER -e $EMAIL_LIST -r $USER -i $CORE -n $TARGET_CORE --skip-check 2>&1 | tee -a $LOG_FILE
		echo "" | tee -a $LOG_FILE
	#List user specific copy jobs that are completed
	elif [[ "${OPERATION^^}" == "LIST" ]]; then
		echo -e -n "${ORANGE}\nList DB copy JOBS\n\t=> dbcopy-client -a $OPERATION -u $ENDPOINT --src_host $SOURCE_SERVER --tgt_host $TARGET_SERVER -e $EMAIL_LIST -r $USER${NC}\n\n"
		dbcopy-client -a $OPERATION -u $ENDPOINT --src_host $SOURCE_SERVER --tgt_host $TARGET_SERVER -e $EMAIL_LIST -r $USER | grep -v -e "Complete"
		exit 0
	else
		echo -e -n "${ORANGE}Did not recognise input params. Please double check input params and retry!${NC}\n\n"
		exit 1
	fi

done < $INPUT_CORES

exit 0
