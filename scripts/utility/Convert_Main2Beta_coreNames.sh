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

## A short script to programmatically transform main side CoreDBs to their appropriate CoreDB name for moving to Beta
## It takes into account the appropriate assembly accession, annotation source and prefix

MAIN_CORES_IN=$1
SCHEMA_FOR_BETA=$2
CMD_HOST=$3
OPTIONAL_CORE_PREFIX=$4
if [[ -z $MAIN_CORES_IN ]] || [[ -z $SCHEMA_FOR_BETA ]] ||  [[ -z $CMD_HOST ]]; then
echo "Usage Convert_Main2Beta_coreNames.sh <list of Ensembl Main cores DBs> <Beta schema required> <Staging host> <*Optional_beta_core_prefix>"
exit 1
fi
OUTPUT_TRANSFORM="Beta_Cores.out.tsv"
if [[ -f $OUTPUT_TRANSFORM ]]; then
    rm ./$OUTPUT_TRANSFORM
fi
echo -e -n "#EnsemblMain_CoreDB\tEnsembl BetaCoreDB [schema version:${SCHEMA_FOR_BETA}]\n" > $OUTPUT_TRANSFORM

if [[ -f main_to_beta_convert.log ]]; then
    rm ./main_to_beta_convert.log
fi
if [[ -n $OPTIONAL_CORE_PREFIX ]]; then
    OPTIONAL_CORE_PREFIX="${4}_"
    echo "Prepending optional core prefix '${OPTIONAL_CORE_PREFIX}'" | tee main_to_beta_convert.log
    sleep 2
fi

function convert_main_to_beta(){

        local INPUT_CORES=$1
        local SUFFIX=$2
        local STAGING_HOST=$3

    while read CORE
        do
            SP_PROD_NAME=`$STAGING_HOST -D $CORE -Ne "SELECT meta_value FROM meta WHERE meta_key = 'species.production_name';"`
            ASSEMBLY_ACC=`$STAGING_HOST -D $CORE -Ne "SELECT meta_value FROM meta WHERE meta_key = 'assembly.accession';"`
            ALT_ASSEMBLY_ACC=`$STAGING_HOST -D $CORE -Ne "SELECT meta_value FROM meta WHERE meta_key = 'assembly.alt_accession';"`

            if [ $SUFFIX == "gb_core" ]; then
                SOURCE="Genbank";
                SHORT_SUFFIX='gb'
                MSG="$SOURCE suffix($SHORT_SUFFIX)"
            elif [ $SUFFIX == "rs_core" ]; then
                SOURCE="RefSeq";
                SHORT_SUFFIX='rs'
                MSG="$SOURCE suffix($SHORT_SUFFIX)"
            elif [ $SUFFIX == "cm_core" ]; then
                SOURCE="Community";
                SHORT_SUFFIX='cm'
                MSG="$SOURCE suffix($SHORT_SUFFIX)"
            elif [ $SUFFIX == "fb_core" ]; then
                SOURCE="FlyBase";
                SHORT_SUFFIX='fb'
                MSG="$SOURCE suffix($SHORT_SUFFIX)"
            elif [ $SUFFIX == "wb_core" ]; then
                SOURCE="WormBase";
                SHORT_SUFFIX='wb'
                MSG="$SOURCE suffix($SHORT_SUFFIX)"
            elif [ $SUFFIX == "vb_core" ]; then
                SOURCE="VectorBase";
                SHORT_SUFFIX='vb'
                MSG="$SOURCE suffix($SHORT_SUFFIX)"
            else
                SOURCE="Non-anno source suffix"
                SHORT_SUFFIX=''
                MSG="$SOURCE (na)"
            fi

            echo -e -n "\t ---- Converting $MSG cores ----\n"

            NEW_SP_PROD_NAME=`echo $SP_PROD_NAME | sed -E s/_gc[af].+//g`

            if [[ -n $ALT_ASSEMBLY_ACC ]]; then
                LC_ACCESSION=`echo $ALT_ASSEMBLY_ACC | sed 's/GCF_/gcf/' | sed 's/\./v/'`
                ACCESSION_VERSION=`echo $ASSEMBLY_ACC | sed -E 's/GCA_[0-9]+//g' | sed 's/\.//'`
                ALT_ACCESSION_VERSION=`echo $ALT_ASSEMBLY_ACC | sed -E 's/GCF_[0-9]+//g' | sed 's/\.//'`
                if [[ $ACCESSION_VERSION != $ALT_ACCESSION_VERSION ]]; then
                    echo -e -n "\t*** Version difference ($CORE) GCA vs GCF accession: GCA(v$ACCESSION_VERSION) vs GCF(v$ALT_ACCESSION_VERSION)\n" \
                    | tee -a main_to_beta_convert.log
                fi
            else
                LC_ACCESSION=`echo $ASSEMBLY_ACC | sed 's/GCA_/gca/' | sed 's/\./v/'`
            fi

            NEW_BETA_CORE_NAME="${OPTIONAL_CORE_PREFIX}${NEW_SP_PROD_NAME}_${LC_ACCESSION}${SHORT_SUFFIX}_core_${SCHEMA_FOR_BETA}_1"
            echo -e -n "$CORE --> $NEW_BETA_CORE_NAME\n" | tee -a main_to_beta_convert.log
            echo -e -n "$CORE\t$NEW_BETA_CORE_NAME\n" >> $OUTPUT_TRANSFORM
        done < $INPUT_CORES
}

if [[ -f temp_combined_suffix.list.tmp ]]; then
    rm ./temp_combined_suffix.list.tmp
fi

for SUFFIX in gb_core rs_core vb_core cm_core wb_core fb_core
do
    TEMP_CORES_LIST="temp_cores_${SUFFIX}.list.tmp"
    grep -e $SUFFIX $MAIN_CORES_IN > $TEMP_CORES_LIST
    SHORT_SUFFIX=`echo $SUFFIX | sed 's/_core//'`;
    cat $TEMP_CORES_LIST >> temp_combined_suffix.list.tmp
    # Convert all cores with source suffix
    convert_main_to_beta $TEMP_CORES_LIST $SUFFIX $CMD_HOST
done

# Convert all cores without any additional source suffix
grep -v -f temp_combined_suffix.list.tmp $MAIN_CORES_IN > non-suffix-cores.tmp
convert_main_to_beta non-suffix-cores.tmp none $CMD_HOST

cat ./temp_combined_suffix.list.tmp | sort > MainRelease_Cores.WithSuffix.txt
cat ./non-suffix-cores.tmp | sort > MainRelease_Cores.NoSuffix.txt
rm ./*.tmp
echo "** See main output TSV -> Beta_Cores.out.tsv **"