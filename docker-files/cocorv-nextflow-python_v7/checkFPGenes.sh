main() {
  CoCoRVFolder=$1
  cocorvOut=$2
  topK=$3
  casecontrol=$4
  build=$5
  sampleList=$6
  addVEP=$7

  vcfSuffix=".biallelic.leftnorm.ABCheck.vcf.gz"
  vcfAnnoSuffix=".annotated.vcf.gz"
  outputFile="top${topK}.${cocorvOut}.${casecontrol}.variants.tsv"
  fullGenotype=T
  
  if [[ ${addVEP} == "T" ]]; then
    annotations="SYMBOL,FILTER,Consequence,HGVSp,REVEL,am_class,am_pathogenicity,SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL,LoF,LoF_filter,LoF_flags,LoF_info,CADD_PHRED,CADD_RAW"
  else  
    annotations="Gene.refGene,FILTER,Func.refGene,ExonicFunc.refGene,AAChange.refGene,REVEL"
  fi

  bash ${CoCoRVFolder}/utilities/postCheckCoCoRV_docker.sh ${cocorvOut} ${topK} \
    ${vcfAnnoSuffix} ${outputFile} ${casecontrol} ${fullGenotype} ${sampleList} ${build} ${annotations} ${vcfSuffix}
}

main "$@"
