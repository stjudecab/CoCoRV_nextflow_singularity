main() {
  CoCoRVFolder=$1
  cocorvOut=$2
  topK=$3
  casecontrol=$4
  cocorvOutputRoot=$5
  build=$6
  sampleList=$7

  variantFolder="${cocorvOutputRoot}/CoCoRV/byChr/"
  vcfPrefix="${cocorvOutputRoot}/vcf_vqsr_normalizedQC/"
  vcfSuffix=".biallelic.leftnorm.ABCheck.vcf.gz"
  vcfAnnoPrefix="${cocorvOutputRoot}/annotation/"
  vcfAnnoSuffix=".annotated.vcf.gz"
  outputFile="${cocorvOut}.${casecontrol}.variants.tsv"
  fullGenotype=T
  annotations="Gene.refGene,FILTER,Func.refGene,ExonicFunc.refGene,AAChange.refGene,REVEL"

  bash ${CoCoRVFolder}/utilities/postCheckCoCoRV_docker.sh ${cocorvOut} ${topK} ${variantFolder} \
    ${vcfAnnoPrefix} ${vcfAnnoSuffix} ${outputFile} ${casecontrol} ${fullGenotype} ${sampleList} ${build} ${annotations} ${vcfPrefix} ${vcfSuffix}
}

main "$@"
