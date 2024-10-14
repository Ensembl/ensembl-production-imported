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

#### A script used to process a bulk Handover of genome cores to production.
# See documentation page: https://github.com/Ensembl/ensembl-prodinf-tools/blob/main/docs/bulk_handover_database.rst

#Input params
HOST=$1
INPUT=$2
MAIN_OR_BETA=${3^}
DIVISION=$4
CONTEXT=$5
CWD=`readlink -f $PWD`

# Check required params first
if [[ -z $HOST ]] || [[ -z $INPUT ]] || [[ -z $MAIN_OR_BETA ]] || [[ -z $DIVISION ]] || [[ -z $CONTEXT ]]; then

	echo "Usage: Bulk_CoreDB_Handover_Beta <Core DB host> <INPUT_CORES_LIST> <main or beta ?> <Division> <'Short HO context description'>"
	echo "E.g. Bulk_CoreDB_Handover_Beta me1 Cores.list.txt main metazoa 'Metazoa handover first run'"
	exit 0
fi

INPUT_CORES_LISTED=`readlink -f $INPUT`

#Ensure core list file isn't empty
if [[ ! -s $INPUT_CORES_LISTED ]]; then
	echo "Input list file of core DB names is empty! Exiting"
	exit 1
fi

# Check an appropriate division was supplied by user
if [[ "$DIVISION" != "metazoa" ]] && [[ "$DIVISION" != "plants" ]] && [[ "$DIVISION" != "microbes" ]] && [[ "$DIVISION" != "vertebrates" ]] && [[ "$DIVISION" != "rapid" ]]; then
	echo "Division supplied ($DIVISION) must be defined as: [ metazoa | plants | microbes | vertebrates | rapid ]"
	exit 1
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

# Find out if we want to do a test run:
read -p "Set run as a TEST dry run ? <yes|no> " DRY_RUN
if [[ "${DRY_RUN^^}" == "YES" ]]; then
	echo -e -n "${ORANGE}SAFE RUN: Handover dry test run.${NC}\n\n"
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
	git clone -b main --depth 1 git@github.com:Ensembl/ensembl-production-imported-private.git $$CWD/ensembl-production-imported-private
	source $CWD/ensembl-production-imported-private/lib/bash/HANDOVER.env
	URL_SOURCE="$CWD/ensembl-production-imported-private/lib/bash/HANDOVER.env"
fi

# Check for endpoint URL correctly defined and set for handover-copy purposes
if [ $HANDOVER_BASE_URL ]; then
	echo -e -n "${GREEN}Detected handover-client endpoint URL variable -> $URL_SOURCE${NC}\n"
	ENDPOINT_BETA="${HANDOVER_BASE_URL}/rapid/handovers/"
	ENDPOINT_MAIN="${HANDOVER_BASE_URL}${DIVISION}/handovers/"
else
	echo -e -n "${RED}Unable to detect Endpoint URL variable for handover-client.\n Required repo: 'ensembl-production-imported-private'${NC}\n\n"
fi

# Check for the right modenv is located and activated
which handover-client
if [ $? -eq 0 ]; then
	echo -e -n "${GREEN}ensembl_mod_env module environment detected 'ensembl/production-tools' handover-client found! Proceeding...${NC}\n\n"
else
	echo -e -n "${RED}Unable to locate 'handover-client' (production-tools). Ensure you have created/loaded via 'module load ensembl/production-tools'.\n"
	echo -e -n "modenv_create ensembl/production-tools production-tools && module load ensembl/production-tools'${NC}\n"
	echo -e -n "${ORANGE}See --> https://github.com/Ensembl/ensembl-mod-env/wiki/\n\n"
	exit 1
fi

# Database host expansion:
DATABASE_SERVER=$($HOST details url)


#Setting End point for RR vs MainSite handover:
if [[ $MAIN_OR_BETA == "Main" ]]; then
	echo -e -n "DB Core handover [MainSite] endpoint --> \"$ENDPOINT_MAIN\"\n"
	SUB_ENDPOINT=$ENDPOINT_MAIN

	#### Find out what Ensembl and EnsemblGenomes release version are needed:
	read -p "Enter now Ensembl release version. e.g. 110: --> " E_RELEASE
	read -p "Enter now EnsemblGenomes release version. e.g. 57: --> " EG_RELEASE
	CONTEXT_DESC="HO Ensembl ${DIVISION^} '${CONTEXT}' (e:${E_RELEASE}, EG:${EG_RELEASE}:"
	RELEASE_DESC="e${E_RELEASE}_EG${EG_RELEASE}"

elif [[ $MAIN_OR_BETA == "Beta" ]]; then
	echo -e -n "DB Core handover [Beta] endpoint -> \"$ENDPOINT_BETA\"\n"
	SUB_ENDPOINT=$ENDPOINT_BETA
	read -p "Enter now Beta/MVP version. e.g. 45 : " BETA_VERSION
	CONTEXT_DESC="Ensembl ${DIVISION^} Beta/MVP HO (Beta/MVP${BETA_VERSION}:"
	RELEASE_DESC="Beta_MVP:${BETA_VERSION}"
else
	echo -e -n "${RED}Can't decide which endpoint to submit core handover. Opts are: 'Main' or 'Beta' !\nExiting${NC}\n"
	exit 1
fi

## Start HO across set of input core(s)
for DB in $(cat $INPUT_CORES_LISTED); do
	SPECIES=$(echo $DB | perl -pe "s/_core_.*//");
	DESCRIPTION="${CONTEXT_DESC} $SPECIES)"
	SOURCE_URL=$($HOST details url $DB)
	if [[ "${DRY_RUN^^}" == "YES" ]]; then
		echo -e -n "${ORANGE}\nDry-run handover test:${NC}\n\n"
		echo -e -n "${ORANGE}handover-client --action submit --uri ${SUB_ENDPOINT} --src_uri \"${DATABASE_SERVER}${DB}\" --email \"${USER}@ebi.ac.uk\" --description \"${DESCRIPTION}\"${NC}\n"
		echo "handover-client --action submit --uri ${SUB_ENDPOINT} --src_uri \"${DATABASE_SERVER}${DB}\" --email \"${USER}@ebi.ac.uk\" --description \"${DESCRIPTION}\"" > ${SPECIES}_HandOver-DRYRUN.${RELEASE_DESC}.${DIVISION}.log
	elif [[ "${DRY_RUN^^}" == "NO" ]]; then
		echo -e -n "${GREEN}\nPerforming real handover: $DB${NC}\n\n"
		handover-client \
			--action submit \
			--uri ${SUB_ENDPOINT} \
			--src_uri "${DATABASE_SERVER}${DB}" \
			--email "${USER}@ebi.ac.uk" \
			--description "${DESCRIPTION}" 2>&1 | tee ${SPECIES}_HandOver-REAL.${RELEASE_DESC}.${DIVISION}.log

		if [[ $? == 0 ]]; then
			echo -e -n "${GREEN}Handover done - $SPECIES${NC}\n"
		else
			echo "${RED}handover-client submission appears to have failed: $SPECIES${NC}"
		fi
		sleep 2
	else
		echo -e -n "${ORANGE}Did not recognise 'Dry Run' param. Options are: <yes|no>.${NC}\n\n"
	fi
	echo ""
done

exit 0

