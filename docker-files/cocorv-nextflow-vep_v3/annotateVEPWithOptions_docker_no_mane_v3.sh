set -euo pipefail

main() {
# usage: bash <script> ${vcfFile} ${ASSEMBLY} ${outputPrefix} ${REF} ${LOFTEEDIR} ${LOFTEEDATADIR} ${CADDSNVS} ${CADDINDELS} ${SPLICEAISNVS} ${SPLICEAIINDELS} ${AlphaMissenseDIR} ${REVELDIR} ${nThreads} ${annotations}


### HANDLE INPUTS
vcfFile=$1
ASSEMBLY=$2
outputPrefix=$3
REF=$4
CADDSNVS=$5
CADDINDELS=$6
SPLICEAISNVS=$7
SPLICEAIINDELS=$8
AlphaMissenseDIR=$9
REVELDIR=${10}
nThreads=${11}
annotations=${12}
VEPCACHE=${13}

### SET VARIABLES
if [ $ASSEMBLY == 'GRCh38' ]; then
  LOFTEE="LoF,loftee_path:${VEPCACHE}/Plugins,human_ancestor_fa:/${VEPCACHE}/Plugins/data/human_ancestor.fa.gz,conservation_file:/${VEPCACHE}/Plugins/data/loftee.sql,gerp_bigwig:/${VEPCACHE}/Plugins/data/gerp_conservation_scores.homo_sapiens.GRCh38.bw --dir_plugins ${VEPCACHE}/Plugins"
  AlphaMissense="AlphaMissense,file=${AlphaMissenseDIR}/hg38/AlphaMissense_hg38.tsv.gz"
  REVEL="REVEL,file=${REVELDIR}/new_tabbed_revel_grch38.tsv.gz,no_match=1"
elif [ $ASSEMBLY == 'GRCh37' ]; then
    LOFTEE="LoF,loftee_path:${VEPCACHE}/Plugins,human_ancestor_fa:${VEPCACHE}/Plugins/data/human_ancestor.fa.gz,conservation_file:${VEPCACHE}/Plugins/data/phylocsf_gerp.sql,gerp_file:${VEPCACHE}/Plugins/data/GERP_scores.final.sorted.txt.gz"
  AlphaMissense="AlphaMissense,file=${AlphaMissenseDIR}/hg19/AlphaMissense_hg19.tsv.gz"
  REVEL="REVEL,file=${REVELDIR}/new_tabbed_revel.tsv.gz,no_match=1"
else
    echo "GENOME VERSION PROVIDED IS INVALID!"
    exit 1;
fi

CADD="CADD,${CADDSNVS},${CADDINDELS}"
SPLICEAI="SpliceAI,snv=${SPLICEAISNVS},indel=${SPLICEAIINDELS}"

#perlThreadForLoftee="/research_jude/rgs01_jude/groups/cab/projects/Control/common/reference/VEP_v103/ensembl-vep/cache_hg38/cpanm/lib/perl5/x86_64-linux-thread"
export PERL5LIB=${PERL5LIB}:${VEPCACHE}/Plugins
echo ${PERL5LIB}

annotationList=(${annotations//,/ })
pluginString=""
columnString=""
for annotation in "${annotationList[@]}"; do
  plugin=""
  column=""
  if [[ ${annotation} == "LOFTEE" ]]; then
    plugin="--plugin ${LOFTEE}"
    column="LoF,LoF_filter,LoF_flags,LoF_info"
  elif [[ ${annotation} == "CADD" ]]; then
    plugin="--plugin ${CADD}"
    column="CADD_PHRED,CADD_RAW"
  elif [[ ${annotation} == "SPLICEAI" ]]; then
    plugin="--plugin ${SPLICEAI}"
    column="SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL"
  elif [[ ${annotation} == "AM" ]]; then
    plugin="--plugin ${AlphaMissense}"
    column="am_class,am_pathogenicity"
  elif [[ ${annotation} == "REVEL" ]]; then
    plugin="--plugin ${REVEL}"
    column="REVEL"
  else
    echo "not supported annotation: ${annotation}"
    exit 1
  fi

  pluginString="${pluginString} ${plugin}"
  if [[ ${columnString} == "" ]]; then
    columnString="${column}"
  else 
    columnString="${columnString},${column}"
  fi
done

### RUN VEP
vep --cache --dir ${VEPCACHE} \
  --offline --vcf --assembly ${ASSEMBLY} \
  ${pluginString} \
  --force_overwrite \
  --compress_output bgzip \
  -i ${vcfFile} \
  -o ${outputPrefix}.tmp.vcf.gz \
  --fork ${nThreads} \
  --pick \
  --variant_class \
  --fasta ${REF} 2> ${outputPrefix}_vep.log


## IMPROVE VEP ANNOTATION BY SPLITTING
bcftools +split-vep \
  -c SYMBOL,Consequence,${columnString} \
  -Oz -o ${outputPrefix}.tmp.vep.vcf.gz ${outputPrefix}.tmp.vcf.gz
tabix -p vcf ${outputPrefix}.tmp.vep.vcf.gz

### APPLY CLEANED ANNOTATIONS
bcftools annotate -a ${outputPrefix}.tmp.vep.vcf.gz \
  --collapse none \
  -c CSQ,SYMBOL,Consequence,${columnString} \
  -Oz -o ${outputPrefix}.vcf.gz ${vcfFile}
tabix -p vcf ${outputPrefix}.vcf.gz

### CLEAN UP
rm ${outputPrefix}.tmp.vcf.gz* ${outputPrefix}.tmp.vep.vcf.gz*

}

main "$@"
