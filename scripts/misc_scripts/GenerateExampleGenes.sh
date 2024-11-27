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

#!/usr/bin/env bash

INPUT_CORES=$1
MODE=${2^^}
PATCH_OUTFILE="${PWD}/$3"
HOST=$4
OUTPUT_TSV_SUMMARY="${PWD}/New_EG_Genes_Selected.tsv"
COMPARA_DB=""

if [[ -z $INPUT_CORES ]] || [[ -z $MODE ]] || [[ -z $PATCH_OUTFILE ]] || [[ -z $HOST ]] || [[ ! -n "$COMPARA_DB" ]]; then
    echo "Usage: sh GenerateExampleGenes.sh <INPUT FILE> <auto|manual> <PATCH_OUTFILE_NAME> <CORE HOST>"
    echo "Input file format"
    echo "Manual mode: core_db_name<'\n'>"
    echo "Auto mode: core_db_name<'\t'>manually_selected_stable_id<'\n'>"
    echo "Make sure to alter script and add the appropriate compara DB!! Line: 22."
    exit 0;
fi

META_OUTDIR=$PWD/GENE_META
mkdir -p $META_OUTDIR

# Some different filters to search for a sub sample of core stable_ids
BIOTYPE="protein_coding" # Not, this really can't be changed as e.g genes are related to peptide compara
STABLE_ID_SQL_FILTER="and description NOT LIKE '%uncharacterized%' and description IS NOT NULL " # SUGGEST TO USE THIS FILTER ON GENE RANDOM SAMPLE 
# STABLE_ID_SQL_FILTER="and description NOT LIKE '%uncharacterized%' "
# STABLE_ID_SQL_FILTER="and description IS NOT NULL "

# Get meta data on core + selected stable_id
gather_gene_meta () {
    local PROD=$1
    local INPUT_STABLE_ID=$2

    echo "Gathering meta data on stable_id: $INPUT_STABLE_ID"
    $HOST -D $CORE -Ne "select stable_id, seq_region_id, seq_region_start, seq_region_end, description from gene where stable_id=\"$INPUT_STABLE_ID\";" \
        > $META_OUTDIR/gene_${PROD}_meta_${INPUT_STABLE_ID}.tmp
    $HOST -D $CORE -Ne "select t.stable_id, t.seq_region_id, t.seq_region_start, t.seq_region_end from transcript t join gene g where t.gene_id = g.gene_id and g.stable_id = \"$INPUT_STABLE_ID\" and g.canonical_transcript_id = t.transcript_id;" \
        > $META_OUTDIR/transcript_${PROD}_meta_${INPUT_STABLE_ID}.tmp
    $HOST -D $CORE -Ne "select s.name from seq_region s join gene g where g.seq_region_id = s.seq_region_id and g.stable_id=\"$INPUT_STABLE_ID\";" \
        > $META_OUTDIR/seq_region_${PROD}_meta_${INPUT_STABLE_ID}.tmp
    $HOST -D $CORE -Ne "select x.display_label from xref x join gene g where x.xref_id = g.display_xref_id and g.stable_id=\"$INPUT_STABLE_ID\";" \
        > $META_OUTDIR/xref_display_${PROD}_meta_${INPUT_STABLE_ID}.tmp

    return
}

if [[ -e $PWD/$OUTPUT_TSV_SUMMARY ]]; then echo "Removing old summary TSV"; rm $OUTPUT_TSV_SUMMARY; fi
echo -e -n "DB Name\tNew Example Gene\tGeneDescription\n" > $OUTPUT_TSV_SUMMARY

if [[ -e $PATCH_OUTFILE ]]; then echo "Removing old PATCH SQL"; rm $PATCH_OUTFILE; fi

while read LINE
do
    CORE=`echo "$LINE" | cut -d $'\t' -f1`
    PROD_NAME=`echo $CORE | sed -E 's/_core.+$//g'`

    if [[ "$MODE" == "AUTO" ]]; then
        echo -e -n "\nAuto patching : $CORE\n"
        $HOST -D $CORE -Ne  "select stable_id from gene where biotype = \"$BIOTYPE\" ${STABLE_ID_SQL_FILTER} order by RAND() limit 25;" \
            | tr "\n" " " >> $META_OUTDIR/rand_genes_${CORE}.txt
        
        NEW_GENE=""
        
        for NEW_GENE_CHECK in $(tail -n 1 $META_OUTDIR/rand_genes_${CORE}.txt);
            do
            # Now check if compara DB contains homologies on randomly select gene stable_id
            CHECK_COMPARA_HOMOLOGY=`$HOST -D $COMPARA_DB -Ne "SELECT count(*) FROM gene_member JOIN genome_db USING(genome_db_id) JOIN gene_member_hom_stats USING(gene_member_id) WHERE genome_db.name = \"${PROD_NAME}\" AND gene_trees > 0 and stable_id = \"${NEW_GENE_CHECK}\";"`

            if [[ $CHECK_COMPARA_HOMOLOGY -ge 1 ]]; then

                echo "Homologies for $NEW_GENE_CHECK located! [$CORE]"

                CHECK_DESCRIPTION=`$HOST -D $CORE -Ne "SELECT COUNT(*) FROM gene WHERE stable_id = \"$NEW_GENE_CHECK\" $STABLE_ID_SQL_FILTER;"`

                if [[ $CHECK_DESCRIPTION -eq 1 ]]; then
                    NEW_GENE=$NEW_GENE_CHECK
                    # Get meta on stable id:
                    gather_gene_meta $PROD_NAME $NEW_GENE_CHECK
                    break
                fi
            else
                echo "Failed homology check: $CORE -> $NEW_GENE_CHECK" | tee -a Missing_homologies.txt
            fi
        done
    elif [[ "$MODE" == "MANUAL" ]]; then

        echo -e -n "\nManual patching : $CORE\n"
        # Quickly take the gene stable_id from the input file 
        NEW_GENE=`echo "$LINE" | cut -d $'\t' -f2`
        gather_gene_meta $CORE $NEW_GENE

    fi

    GENE_METAFILE="$META_OUTDIR/gene_${PROD_NAME}_meta_${NEW_GENE}.tmp"
    SEQR_METAFILE="$META_OUTDIR/seq_region_${PROD_NAME}_meta_${NEW_GENE}.tmp"
    TRANSC_METAFILE="$META_OUTDIR/transcript_${PROD_NAME}_meta_${NEW_GENE}.tmp"
    XREF_METAFILE="$META_OUTDIR/xref_display_${PROD_NAME}_meta_${NEW_GENE}.tmp"

    GENE_DESCRIPTION=`cut -f5 "$GENE_METAFILE"`
    echo -e -n "$CORE\t$NEW_GENE\t$GENE_DESCRIPTION\n" >> $OUTPUT_TSV_SUMMARY

    # # gene info:
    GENE_STABLE_ID=$NEW_GENE
    SEQR_ID=`cut -f2 $GENE_METAFILE`
    SEQR_NAME=`cut -f1 $SEQR_METAFILE`
    GENE_SEQR_START=`cut -f3 $GENE_METAFILE`
    GENE_SEQR_END=`cut -f4 $GENE_METAFILE`

    # gene info:
    TRANSCRIPT_STABLE_ID=`cut -f1 $TRANSC_METAFILE`
    TRANSCRIPT_TEXT=$TRANSCRIPT_STABLE_ID
    TRANSCRIPT_SEQR_START=`cut -f3 $TRANSC_METAFILE`
    TRANSCRIPT_SEQR_END=`cut -f4 $TRANSC_METAFILE`

    #XREF info:
    XREF_DISPLAY=`cut -f1 $XREF_METAFILE`

    if [[ "$XREF_DISPLAY" == "" ]]; then
        XREF_DISPLAY=$GENE_STABLE_ID
    fi

    # Write new meta data to Patch sql outfile
    echo "USE $CORE;" >> $PATCH_OUTFILE
    echo "UPDATE meta SET meta_value = \"$NEW_GENE\" where meta_key = \"sample.gene_param\";" >> $PATCH_OUTFILE
    echo "UPDATE meta SET meta_value = \"$XREF_DISPLAY\" where meta_key = \"sample.gene_text\";" >> $PATCH_OUTFILE
    echo "UPDATE meta SET meta_value = \"${SEQR_NAME}:${GENE_SEQR_START}-${GENE_SEQR_END}\" where meta_key = \"sample.location_param\";" >> $PATCH_OUTFILE
    echo "UPDATE meta SET meta_value = \"$TRANSCRIPT_STABLE_ID\" where meta_key = \"sample.transcript_param\";" >> $PATCH_OUTFILE
    echo "UPDATE meta SET meta_value = \"${SEQR_NAME}:${TRANSCRIPT_SEQR_START}-${TRANSCRIPT_SEQR_END}\" where meta_key = \"sample.location_text\";" >> $PATCH_OUTFILE
    echo -e -n "UPDATE meta SET meta_value = \"$TRANSCRIPT_TEXT\" where meta_key = \"sample.transcript_text\";\n\n" >> $PATCH_OUTFILE

done < $INPUT_CORES
