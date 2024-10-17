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

############################################################################
# This script function is to gather and preserve assembly snapshots of a given Ensembl Release. 
# It operates in three main stages:
# 1) - Gather all Core databases on a specific staging host
# 2) - Parse all core DBs by Ensembl divisions
# 3) - Iterate over one specific division defined by the user. 
# 4) - Parse key assembly/annotation related meta table keys and create 
#      a combined snapshot of all DBs meta into a TSV outputfile.
# 5) - Gather related taxonomy information per DB and output this to an additional TSV output file.
# 6) - Outputs some summary info on databases processed and taxonomy info.
# 7) - Compare and display info on current release to a previous release snapshot for DB changes/additions 
#      (If available i.e. already generated)

# Typical usage: Post Ensembl Core DB Handover; run this script on desired staging host (e.g. st3b) 
# or pre-staging host (e.g. st3).
# If available, already generated TSVs can be compared between Ensembl releases. 
############################################################################

## Sript version: v1.0
## Created: Lahcen Campbell (Ensembl Metazoa) - [lcampbell@ebi.ac.uk]

DIVISION=$1
RELEASE_HOST=$2
RELEASE=$3
PREVIOUS_HOST=$4
PREVIOUS_RELEASE=$5
CWD=`readlink -f $PWD`
DATE=`date | awk {'print $1,$2,$3,$6'} | sed 's/ /_/g'`

# ## Taxonomy information parsing.
TAXONOMY_SCRIPT="${ENSEMBL_ROOT_DIR}/ensembl-production-imported/scripts/assembly_release_tracking/GetTaxonomy.pl"

# Check for minimum information required to process a given release.
if [[ -z $RELEASE_HOST ]] || [[ -z $RELEASE ]] || [[ -z $DIVISION ]]; then
	echo "Usage: Ensembl_MainRelease_Tracker.sh <Division: metazoa|plants|fungi|vertebrates|protists> <staging host: st1, st3, st4, st1b, st3b, st4b> <Ensembl Version> Optional:<Past staging Host> <Past version>" #e.g. st3 108 st3b 107
	echo -e -n "E.g: Ensembl_MainRelease_Tracker metazoa st3 114 st3b 113\n"
	echo -e -n "     Ensembl_MainRelease_Tracker vertebrates st1 113\n\n"
	echo -e -n "Note: Your working directory for this pipeline should contain a clone of 'ensembl-production-imported'!\n"
	exit 1
fi

##Vars for STDOUT colour formatting
read -p "Do you want colourised STDOUT printing enabled ? <yes|no> " COLOUR_OPT
if [[ "$COLOUR_OPT" == "YES" ]] || [[ "$COLOUR_OPT" == "yes" ]]; then
	GREEN='\033[0;32m'
	RED='\033[0;31m'
	ORANGE='\033[0;33m'
	PURPLE='\033[0;35m'
	NC='\033[0m' # <- No Colour
else
	GREEN=''
	RED=''
	ORANGE=''
	PURPLE=''
	NC=''
fi

#Expand release host from codon variable to full host:port for taxonomy registry:
RELEASE_PORT=`echo $($RELEASE_HOST details script) | cut -d ' ' -f4`
TAXON_RELEASE_HOST=`echo $($RELEASE_HOST details script) | cut -d ' ' -f2`

echo -e -n "\n${PURPLE}Creating database snapshot of host:$RELEASE_HOST (port:$RELEASE_PORT). Ensembl release version: $RELEASE\n\n"

# If previous release information not given, generate default settings.
if [[ -z $PREVIOUS_RELEASE ]]; then
	PREVIOUS_RELEASE=$((RELEASE-1))

	if [[ -z $PREVIOUS_HOST ]] && [[ $RELEASE_HOST != +(st1b|st3b|st4b) ]]; then
	PREVIOUS_HOST="${RELEASE_HOST}b"
		echo -e -n "${ORANGE}Previous staging host was not provided !\nBest guess is to use 'Odd' numbered staging host for release snapshot comparison (if available) - host:$PREVIOUS_HOST, version:$PREVIOUS_RELEASE${NC}\n\n"
	elif [[ -z $PREVIOUS_HOST ]] && [[ $RELEASE_HOST == +(st1b|st3b|st4b) ]]; then
		PREVIOUS_HOST=`echo ${RELEASE_HOST} | sed 's/b//'`
		echo -e -n "${ORANGE}Previous staging host was not provided !\nBest guess is to use 'Even' numbered staging host for release snapshot comparison (if available) - host:$PREVIOUS_HOST, version:$PREVIOUS_RELEASE${NC}\n\n"
	fi
else
	echo -e -n "Will compare release snapshot with that of previous host $PREVIOUS_HOST on Ensembl release: $PREVIOUS_RELEASE${NC}\n\n"
fi

# Check an appropriate division was supplied by user
if [[ "$DIVISION" != "metazoa" ]] \
	&& [[ "$DIVISION" != "plants" ]] \
	&& [[ "$DIVISION" != "fungi" ]] \
	&& [[ "$DIVISION" != "vertebrates" ]] \
	&& [[ "$DIVISION" != "protists" ]]; then
	echo -e -n "${RED}Division supplied ($DIVISION) not recognised. Must define as: [ metazoa | plants | fungi | vertebrates | protists ]${NC}\n\n"
	exit 1
fi

#Create folder for current release and gather files from past release
CUR_REL_FOLDER="${CWD}/E${RELEASE}"
PREV_REL_FOLDER="${CWD}/E${PREVIOUS_RELEASE}"
mkdir -p $CUR_REL_FOLDER
cd $CUR_REL_FOLDER

# Main Output files
ALL_CORE_DBS="${RELEASE_HOST}_all_e${RELEASE}_coredbs.txt"
METAZOA_CORES="${RELEASE_HOST}_metazoa_cores_${RELEASE}.txt"
FUNGI_CORES="${RELEASE_HOST}_fungi_cores_${RELEASE}.txt"
PLANTS_CORES="${RELEASE_HOST}_plants_cores_${RELEASE}.txt"
VERTEBRATES_CORES="${RELEASE_HOST}_vertebrates_cores_${RELEASE}.txt"
PROTISTS_CORES="${RELEASE_HOST}_protists_cores_${RELEASE}.txt"

# Final and temp Output files
ASSEMBLY_INFO="${CUR_REL_FOLDER}/${DIVISION^}_sp_Asm_${RELEASE_HOST}_e${RELEASE}.info.tsv"
PAST_ASSEMBLY_INFO="${PREV_REL_FOLDER}/${DIVISION^}_sp_Asm_${PREVIOUS_HOST}_e${PREVIOUS_RELEASE}.info.tsv"
CHANGES_BETWEEN_REL="${DIVISION^}_species_Diff_e${RELEASE}_comparedTo_e${PREVIOUS_RELEASE}.list.txt"

TEMP_SNAPSHOT_FILE="${CUR_REL_FOLDER}/temp_snapshot.tsv"
if [[ -e $TEMP_SNAPSHOT_FILE ]]; then rm $TEMP_SNAPSHOT_FILE; fi

## Stage 1a - get all core databases from staging host
if [[ -s $ALL_CORE_DBS ]]; then

	ALL_CORE_COUNT=$(wc -l < $ALL_CORE_DBS)
	echo -e -n "Located previously generated core list: $ALL_CORE_DBS\n\t${ORANGE}#> Contains $ALL_CORE_COUNT core dbs from e${RELEASE}.${NC}\n"

elif [[ ! -z $RELEASE_HOST ]]; then

	if [[ $DIVISION == "vertebrates" ]]; then
		CORE_REGEX="core_${RELEASE}"
	else
		CORE_REGEX="core_.+${RELEASE}"
	fi

	echo -e -n "Finding set of Ensembl ${DIVISION^} cores in host: $RELEASE_HOST for release e${RELEASE}\n"
	$RELEASE_HOST -e "SHOW DATABASES;" | grep -v "Database" | grep -v -e "collection" | grep -e 'core' | grep -E "$CORE_REGEX" > $ALL_CORE_DBS
	ALL_CORE_COUNT=$(wc -l < $ALL_CORE_DBS)
	echo -e -n "> $ALL_CORE_DBS contains $ALL_CORE_COUNT core dbs in total. <\n\n"

else
	echo "Exiting Stage 1a. Something not quite right !!"
	exit 1
fi

##Stage 1b 
# Set the appropriate target division cores
if [[ "$DIVISION" == "metazoa" ]]; then
	TARGET_DIVISION_CORES=$METAZOA_CORES
elif [[ "$DIVISION" == "plants" ]]; then
	TARGET_DIVISION_CORES=$PLANTS_CORES
elif [[ "$DIVISION" == "fungi" ]]; then
	TARGET_DIVISION_CORES=$FUNGI_CORES
elif [[ "$DIVISION" == "vertebrates" ]]; then
	TARGET_DIVISION_CORES=$VERTEBRATES_CORES
elif [[ "$DIVISION" == "protists" ]]; then
	TARGET_DIVISION_CORES=$PROTISTS_CORES
fi

if [[ -s ${CUR_REL_FOLDER}/$TARGET_DIVISION_CORES ]]; then

	echo -e -n "Skipping stage to locate ** ${DIVISION^} ** cores in $RELEASE_HOST\n\t"
	TARGET_CORE_COUNT=$(wc -l < $TARGET_DIVISION_CORES)
	echo -e -n "${ORANGE}#> Located $TARGET_CORE_COUNT ${DIVISION^} cores from e${RELEASE}${NC}\n\n"

else

	echo "Generating list of Ensembl ${DIVISION^} core dbs:"

	# Assertain which cores dbs belong to which division
	while read CORE
	do
		SQL_DIVISION=`$RELEASE_HOST -D $CORE -e "SELECT meta_value FROM meta WHERE meta_key = 'species.division';" | tail -n 1`
		#METAZOA
		if [[ $SQL_DIVISION == 'EnsemblMetazoa' ]]; then
			echo "Metazoa --> $CORE"
			echo "$CORE" >> $METAZOA_CORES
		#FUNGI
		elif [[ $SQL_DIVISION == 'EnsemblFungi' ]]; then
			echo "Fungal --> $CORE"
			echo "$CORE" >> $FUNGI_CORES
		#PLANTS
		elif [[ $SQL_DIVISION == 'EnsemblPlants' ]]; then
			echo "Plant --> $CORE"
			echo "$CORE" >> $PLANTS_CORES
		#PROTISTS
		elif [[ $SQL_DIVISION == 'EnsemblVertebrates' ]]; then
			echo "Vertebrate --> $CORE"
			echo "$CORE" >> $VERTEBRATES_CORES
		#VERTEBRATES
		elif [[ $SQL_DIVISION == 'EnsemblProtists' ]]; then
			echo "Protist --> $CORE"
			echo "$CORE" >> $PROTISTS_CORES
		fi
	done < $ALL_CORE_DBS

	## critical for proper function: check staging host for vert or non-vert specific staging:
	if [[ $RELEASE_HOST == +(st3|st4|st3b|st4b) ]]; then
		METAZOA_CORE_COUNT=$(wc -l < $METAZOA_CORES)
		FUNGI_CORE_COUNT=$(wc -l < $FUNGI_CORES)
		PROTIST_CORE_COUNT=$(wc -l < $PROTISTS_CORES)
		PLANT_CORE_COUNT=$(wc -l < $PLANTS_CORES)
		VERT_CORE_COUNT=0
		DISPLAY_HOST="nonverts"
	else
		VERT_CORE_COUNT=$(wc -l < $VERTEBRATES_CORES)
		METAZOA_CORE_COUNT=0
		FUNGI_CORE_COUNT=0
		PROTIST_CORE_COUNT=0
		PLANT_CORE_COUNT=0
		DISPLAY_HOST="vertebrate"
	fi
	
	# Now depending on staging host, display number of cores per each division
	if [[ "$DISPLAY_HOST" == "nonverts" ]]; then
		echo -e -n "${GREEN}\nLocated:\n\t> $METAZOA_CORE_COUNT EnsemblMetazoa cores.\n\t> $FUNGI_CORE_COUNT EnsemblFungi cores."
		echo -e -n "\n\t> $PROTIST_CORE_COUNT EnsemblProtist cores.\n\t> $PLANT_CORE_COUNT EnsemblPlants cores.\n\t **Vertebrate cores not considered stage host != 'st1|st1b'${NC}\n\n"
	elif [[ "$DISPLAY_HOST" == "vertebrate" ]]; then
		echo -e -n "${GREEN}\nLocated:\n\t> $VERT_CORE_COUNT EnsemblVertebrates cores.\n\t **Non-vertebrate cores not considered stage host != 'st3|st3b|st4|st4b'${NC}\n\n"
	fi

	# Now set which target databases should be processed based on user param for division
	if [[ "$DIVISION" == "metazoa" ]]; then
		TARGET_DIVISION_CORES=$METAZOA_CORES
		TARGET_CORE_COUNT=$METAZOA_CORE_COUNT
	elif [[ "$DIVISION" == "plants" ]]; then
		TARGET_DIVISION_CORES=$PLANTS_CORES
		TARGET_CORE_COUNT=$PLANT_CORE_COUNT
	elif [[ "$DIVISION" == "fungi" ]]; then
		TARGET_DIVISION_CORES=$FUNGI_CORES
		TARGET_CORE_COUNT=$FUNGI_CORE_COUNT
	elif [[ "$DIVISION" == "vertebrates" ]]; then
		TARGET_DIVISION_CORES=$VERTEBRATES_CORES
		TARGET_CORE_COUNT=$VERT_CORE_COUNT
	elif [[ "$DIVISION" == "protists" ]]; then
		TARGET_DIVISION_CORES=$PROTISTS_CORES
		TARGET_CORE_COUNT=$PROTIST_CORE_COUNT
	fi
fi



## Stage 2 - Gather core information and print to file
echo -e -n "${PURPLE}Entering Assembly status snapshot stage....${NC}\n\n"
if [[ -s $ASSEMBLY_INFO ]]; then
	echo -e -n "${GREEN}#> Located previously generated ${DIVISION^} e${RELEASE} snapshot TSV: $ASSEMBLY_INFO.\n${ORANGE}Skipping DB meta table query stage.${NC}\n\n"

	PROCESSED=$(wc -l < $ASSEMBLY_INFO)
	REAL_PROCESSED=$(wc -l < $ASSEMBLY_INFO)
	REAL_PROCESSED=$((PROCESSED-1))
fi


## Stage 3 - Iterate over all target division cores, check if the TSV is already completed or parse/generate meta information

# Some params to track parsing vs. de-novo generation of meta table info
PARSE_LOCK="NO"
PARSE_COUNT=0
GEN_COUNT=0

while read DATABASE_NAME; do

	if [[ -e $ASSEMBLY_INFO ]] && [[ $REAL_PROCESSED -eq $TARGET_CORE_COUNT ]]; then

		echo -e -n "${PURPLE}All $TARGET_CORE_COUNT e${RELEASE} ${DIVISION^} cores appear to be processed. Job done already YAY!!${NC}"
		TSV_COMPLETE="ALLDONE"
		break 2
	else
		#If meta info file is already generated don't regenerate
		if [[ -s ${CUR_REL_FOLDER}/${DATABASE_NAME}_meta.info || -s ${CUR_REL_FOLDER}/META_INFO_${RELEASE}/${DATABASE_NAME}_meta.info ]]; then
			
			# Meta '.info' file  ${DATABASE_NAME}_meta.info exists.
			if [[ $PARSE_LOCK == "NO" ]]; then
				echo -e -n "${ORANGE}Parsing existing DB Meta files${NC}\n."
				PARSE_LOCK="YES"
				((PARSE_COUNT++))
			else
				echo -e -n "."
				((PARSE_COUNT++))
			fi

		#Otherwise generate the meta info file via query to host+db	
		else #Generate the meta info file by querying main release HOST
			echo -e -n "\nGenerating meta info file: ${DATABASE_NAME}_meta.info\n"
			$RELEASE_HOST -D $DATABASE_NAME -Ne "SELECT meta_key,meta_value FROM meta WHERE meta_key IN ('annotation.provider_name','assembly.provider_name','assembly.default','assembly.accession','genebuild.version','species.scientific_name','species.taxonomy_id');" > ${DATABASE_NAME}_meta.info
			((GEN_COUNT++))
			PARSE_LOCK="NO"
		fi

		METAFILE=`find $CUR_REL_FOLDER -maxdepth 2 -type f -name "${DATABASE_NAME}_meta.info"`
		ANNO_PROVIDER=`grep -w -e "annotation.provider_name" ${METAFILE} | cut -f2 | tr -d '\n'`
		if [[ $ANNO_PROVIDER == "" ]]; then ANNO_PROVIDER="N/A"; fi
		GENEBUILD_VERSION=`grep -w -e "genebuild.version" ${METAFILE} | cut -f2 | tr -d '\n'`
		if [[ $GENEBUILD_VERSION == "" ]]; then GENEBUILD_VERSION="N/A"; fi
		ASM_PROVIDER=`grep -w -e "assembly.provider_name" ${METAFILE} | cut -f2 | tr -d '\n'`
		ASM_DEFAULT=`grep -w -e "assembly.default" ${METAFILE} | cut -f2 | tr -d '\n'`
		ASM_ACCESSION=`grep -w -e "assembly.accession" ${METAFILE} | cut -f2 | tr -d '\n'`
		if [[ $ASM_ACCESSION == "" ]]; then echo -e -n "\n${RED}${DATABASE_NAME} missing 'assembly.accession' meta_value \n${NC}"; fi
		GENUS_SP_NAME=`grep -w -e "species.scientific_name" ${METAFILE} | cut -f2 | tr -d '\n'`
		TAXON_ID=`grep -w -e "species.taxonomy_id" ${METAFILE} | cut -f2`

		echo -e -n "$GENUS_SP_NAME\t$TAXON_ID\t$ANNO_PROVIDER\t$GENEBUILD_VERSION\t$ASM_PROVIDER\t$ASM_DEFAULT\t$ASM_ACCESSION\t$DATABASE_NAME\n" >> $TEMP_SNAPSHOT_FILE
	fi

done < $TARGET_DIVISION_CORES


## Stage 4: Print out basic processing information
echo -e -n "\nIn total [$PARSE_COUNT] existing meta files were parsed.\n"
echo -e -n "In total [$GEN_COUNT] meta files were generated during this run.\n"

## Finalising content and tmp files unless no new core processing performed
if [[ $TSV_COMPLETE != "ALLDONE" ]]; then

	echo -e -n "#Organism\tTaxon ID\tAnnotation provider\tGenebuild version\tAssembly provider\tAsm default\tAsm acc\tCore database\n" > ${ASSEMBLY_INFO}
	cat $TEMP_SNAPSHOT_FILE | sort -k1 >> ${ASSEMBLY_INFO}

	#Locate any and all meta files and move them to appropriate folder
	mkdir -p ${CUR_REL_FOLDER}/META_INFO_${RELEASE}/
	meta_files=( ${CUR_REL_FOLDER}/*_meta.info )
	for MF in "${meta_files[@]}"; do
		if [[ -e "$MF" ]]; then
			mv $MF ${CUR_REL_FOLDER}/META_INFO_${RELEASE}
		fi
	done
	rm $TEMP_SNAPSHOT_FILE
fi

## Stage 5: Generate Taxonomy level information on all cores recorded in the main snapshot TSV 
echo -e -n "${PURPLE}\n\nNow retreiving organismal taxonomic information ...\n"
echo -e -n "Checking for available 'ncbi_taxonomy_${RELEASE}' DB on $RELEASE_HOST...${NC}\n"
CHECK_TAXON_DB=`$RELEASE_HOST -Ne "SHOW DATABASES LIKE 'ncbi_taxonomy_${RELEASE}';"`
TAXONOMY_INPUT_FILE="all_species_and_taxonid.tmp"
TAXONOMY_OUTPUT_FILE="Taxon_levels_${DIVISION}_e${RELEASE}.tsv"
TEMP_TAXON_TSV="taxonomy_input.tsv.tmp"

if [[ $CHECK_TAXON_DB == "ncbi_taxonomy_${RELEASE}" ]]; then
	echo -e -n "${GREEN}Found ncbi taxonomy DB on $RELEASE_HOST.${NC}\n"
	
	# Create or overwrite temp division specific input list of taxa
	cat $ASSEMBLY_INFO | grep -v -e "^#Organism" | cut -f1,2 > $TAXONOMY_INPUT_FILE

	if [[ ! -e $TAXONOMY_OUTPUT_FILE ]]; then

		echo "No previous taxon TSV found, generating afresh now."
		GENERATE_TAXONOMY="YES"

	else
		SNAPSHOT_PROCESSED=$(wc -l < $ASSEMBLY_INFO)
		REAL_PROCESSED=$((SNAPSHOT_PROCESSED -1))
		TAXON_PROCESSED=$(wc -l < $TAXONOMY_OUTPUT_FILE)
		TAXON_REAL_PROCESSED=$((TAXON_PROCESSED -1))

		if [[ $TAXON_REAL_PROCESSED == $REAL_PROCESSED ]]; then
			echo -e -n "${GREEN}\n#> Located previously generated Taxonomy TSV --> '${CUR_REL_FOLDER}/$TAXONOMY_OUTPUT_FILE'.\n${ORANGE}Skipping taxonomy stage.${NC}\n\n"
		else
			echo "Old or partial taxonomy TSV found, generating afresh now."
			GENERATE_TAXONOMY="YES"
		fi

	fi

	#Generate taxonomy info TSV depending on absence or partial taxonomy output file
	if [[ $GENERATE_TAXONOMY == "YES" ]]; then

		# First check for script:
		if [[ ! -e $TAXONOMY_SCRIPT ]]; then echo -e -n "${RED}Couldn't find './ensembl-production-imported/scripts/assembly_release_tracking/GetTaxonomy.pl' ! in $CWD${NC}\n\n"; exit 1; fi

		echo -e -n "#Phylum\tClass\tSuperorder\tOrder\tSuborder\tFamily\tSubfamily\tSpecies\n" > $TAXONOMY_OUTPUT_FILE
		while read TAXON_LINE; do
			echo -e -n "Name\tTaxon ID\n" > $TEMP_TAXON_TSV
			echo "$TAXON_LINE" >> $TEMP_TAXON_TSV
			perl $TAXONOMY_SCRIPT -ncbi_taxon_host $TAXON_RELEASE_HOST -ncbi_taxon_port $RELEASE_PORT -sp_taxon_tsv $TEMP_TAXON_TSV -ens_release $RELEASE | cut -f 1,2,3,4,5,6,7,8 | tail -n 1 >> $TAXONOMY_OUTPUT_FILE
		done < $TAXONOMY_INPUT_FILE
		rm ./$TAXONOMY_INPUT_FILE ./$TEMP_TAXON_TSV

		echo -e -n "${GREEN}Taxonomy TSV generated !${NC}\n\n"
	fi

else
	echo -e -n "${RED}Can't retrieve taxonomy info! Unable to locate database 'ncbi_taxonomy_${RELEASE}' on mysql host: $TAXON_RELEASE_HOST.\nCheck for available ncbi_taxonomy dbs on host.${NC}\n\n"
fi

# exit

### Stage 6:
## Some basic stats on Genus, species and the count of primary haplotype genomes and the distribution of different gene build methods:
SP_COUNT=`cat $ASSEMBLY_INFO | grep -v -e "^#Organism" | cut -f 1 | sort | uniq | wc -l`
GENUS_COUNT=`cat $ASSEMBLY_INFO | grep -v -e "^#Organism" | awk {'print $1'} | sort | uniq | wc -l`
FAMILY_COUNT=`cat $ASSEMBLY_INFO | grep -v -e "^#Organism" | cut -f 6 | sort | uniq | wc -l`

echo -e -n "${PURPLE}### Breakdown of processing ${DIVISION^} release: [ e${RELEASE} | $DATE ] ###\n"
echo -e -n "${NC}#> No. of unique ${DIVISION^} Families = $FAMILY_COUNT\n"
echo -e -n "#> No. of unique ${DIVISION^} Genus = $GENUS_COUNT\n"
echo -e -n "#> No. of unique ${DIVISION^} Species = $SP_COUNT${PURPLE}\n"

echo -e -n "\n#Unique Phyla [Primary assemblies only]:\n"
echo -e -n "Count | Phylum${NC}\n"
cat $TAXONOMY_OUTPUT_FILE | grep -v -e "^#Phylum" -e "alternate_haplotype" | cut -f 1 | sort | uniq -c | sort -nr

echo -e -n "\n${PURPLE}#Unique Orders [Primary assemblies only]:\n"
echo -e -n "Count | Order:${NC}\n"
cat $TAXONOMY_OUTPUT_FILE | grep -v -e "^#Phylum" -e "alternate_haplotype" | cut -f 4 | sort | uniq -c | sort -nr

echo -e -n "\n${PURPLE}#Annotation providers:\n"
echo -e -n "Count | Institute\n${NC}"
cat $ASSEMBLY_INFO | grep -v -e "^#Organism" -e "alternate_haplotype" | cut -f 3 | sort | uniq -c | sort -nr

### Stage 7: 
# Comparsion of TSV files if past release meta TSV file is found
# Checks for past release to find new changes/additions between releases:
if [[ -e $PAST_ASSEMBLY_INFO ]]; then
	echo -e -n "\n${PURPLE}#> Located Meta information for previous release: $PREVIOUS_RELEASE !\nComparing snapshot differences...\n"

	for TSV in $PAST_ASSEMBLY_INFO $ASSEMBLY_INFO; do
		TEMP_SORT_ASM=`echo $TSV | basename $TSV .info.tsv`
		cat $TSV | awk -F"\t" '{print $1,$7}' | grep -v -e "^#Organism" | sort > ${TEMP_SORT_ASM}.sort.tmp
	done

	echo -e -n "\n${ORANGE}# ${DIVISION^} differences in release e$RELEASE compared to release e$PREVIOUS_RELEASE:\n${PURPLE}vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n${NC}"
	comm -23 ${DIVISION^}_sp_Asm_${RELEASE_HOST}_e${RELEASE}.sort.tmp ${DIVISION^}_sp_Asm_${PREVIOUS_HOST}_e${PREVIOUS_RELEASE}.sort.tmp | tee ${CUR_REL_FOLDER}/${CHANGES_BETWEEN_REL}
	echo -e -n "${PURPLE}^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"
	rm ${CUR_REL_FOLDER}/*.sort.tmp 
	sed -Ei 's/ GCA_/\tGCA_/' ${CUR_REL_FOLDER}/${CHANGES_BETWEEN_REL}

	ADDITION_COUNT=$(wc -l < ${CUR_REL_FOLDER}/${CHANGES_BETWEEN_REL})
	echo -e -n "\n${ORANGE}**** $ADDITION_COUNT core addition(s) | change(s) since previous release ****\n"

else
	echo -e -n "\n!!!${ORANGE} Unable to compare release snapshots. No previous release folder: [ E${PREVIOUS_RELEASE} ] or snapshot TSV found.${NC} !!! \n\n"
fi

echo -e -n "\n${PURPLE}Finished processing Ensembl${DIVISION^} release e${RELEASE} snapshot.\n${GREEN}>> $ASSEMBLY_INFO <<\n>> ${CUR_REL_FOLDER}/$TAXONOMY_OUTPUT_FILE <<\n\n${NC}"
exit 0
