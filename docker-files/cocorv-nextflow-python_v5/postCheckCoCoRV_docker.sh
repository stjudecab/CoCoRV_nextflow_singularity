set -eu

main() {
  chrGeneList=$1
  topK=$2
  variantFolder=$3
  vcfAnnoPrefix=$4
  vcfAnnoSuffix=$5
  outputFile=$6
  casecontrol=$7
  fullGenotype=$8
  sampleList=$9
  build=${10}       # GRCh37 or GRCh38
  annotations=${11}  # comma separated annotations

  vcfGTPrefix=""
  vcfGTSuffix=""
  if [[ $# -ge "13" ]]; then
    vcfGTPrefix=${12}
    vcfGTSuffix=${13}
  fi

  # echo ${annotations}

  # make sure the temporary files do not exist
  if [[ -f ${outputFile}.sh || -f ${outputFile}.gt.tsv || -f ${outputFile}.anno.tsv ]]; then
    echo "The following files will be overwritten:"
    echo "${outputFile}.sh"
    echo "${outputFile}.gt.tsv"
    echo "${outputFile}.anno.tsv."
    echo "Move or remove these files first"
    exit 1
  fi

  sampleOption=""
  if [[ ${sampleList} != "NA" ]]; then
    sampleOption="-S ${sampleList}"
  fi
   
  echo -n "" > ${outputFile}
  variantsTempFile=${outputFile}.temp.variants
  addHeader=T

  # header when full genotype is available
  # Gene.refGene,variant,FILTER,QD,Func.refGene,ExonicFunc.refGene,AAChange.refGene,REVEL,AC,AN,AF_popmax
  # header when full genotype is not available
  # Gene.refGene,variant,Func.refGene,ExonicFunc.refGene,AAChange.refGene,REVEL,AC,nhomalt,AN,controls_AC_nfe,AC_nfe,AF_nfe,AF_popmax
  # 

  annotationsHeader=$(echo ${annotations} | sed 's/,/\\t/g')
  annotationsQuery=$(echo ${annotations} | sed 's/,/\\t\%/g')

  while IFS=$'\t' read -r -a myArray; do
    fileID=${myArray[0]}
    gene=${myArray[1]}
    gene=$(echo ${gene} | sed 's/\r//')  # remove \r from Windows's format
    gene=$(echo ${gene} | sed "s/\\\\/\\\\\\\\/g") # in case there are slashes
    if [[ ${fileID} == "" ]]; then
      continue
    fi
    echo "#### ${fileID} ${gene} ####"
    
    chrString=""
    if [[ ! -f ${variantFolder}/${fileID}_snp_indel_file_${casecontrol}.txt && ! -f ${variantFolder}/${fileID}.*${casecontrol}.group ]]; then
      chrString="chr"
    fi
    variantFile=${variantFolder}/${chrString}${fileID}_snp_indel_file_${casecontrol}.txt
    if [[ ! -f ${variantFile} ]]; then
      variantFile=$(ls ${variantFolder}/${fileID}.*${casecontrol}.group); 
    fi

    variants=$(grep "^${gene}\\s" ${variantFile} | sed "s/^${gene}\\s//")
    echo ${variants}

    if [[ ${variants} == "" ]]; then
      continue
    fi

    # should work for both : or - as the separator
    position=$(echo ${variants} | sed 's/:[ATCG]*:[ATCG]*//g' | sed 's/-[ATCG]*-[ATCG]*//g' | sed 's/-/:/g' | sed 's/ /,/g')
    chr=$(echo ${variants} | sed 's/ .*//' | sed 's/:.*//' | sed 's/-.*//')

    if [[ ${build} == "GRCh38" ]]; then
      position=$(echo ${position} | sed -e 's/^/chr/g' -e 's/,/,chr/g')
    fi
    
    echo ${variants} | sed -e 's/ /\n/g' -e 's/,/\n/g' -e's/:/-/g' > ${variantsTempFile}

    if [[ ${fullGenotype} == "T" ]]; then
      if [[ ${addHeader} == "T" ]]; then
        echo -e "variant\tSAMPLE\tGT\tAD\tQD\tAC\tAN\t${annotationsHeader}" > ${outputFile}
        addHeader=F
      fi
      if [[ ${vcfGTPrefix} == "" ]]; then
        # for full genotype cases with annotations
        cmd="bcftools view ${sampleOption} -r ${position} ${vcfAnnoPrefix}${fileID}${vcfAnnoSuffix} | bcftools query -f'[%CHROM-%POS-%REF-%ALT\t%SAMPLE\t%GT\t%AD\t%"${annotationsQuery}"\n]' -i'AC>0 & GT=\"alt\"' | grep -f ${variantsTempFile} | sort -k4,4 >> ${outputFile}"
        # echo $cmd
        echo $cmd > ${outputFile}.sh
        bash ${outputFile}.sh
      else
        # extract genotypes
        bcftools view ${sampleOption} -r ${position} ${vcfGTPrefix}${fileID}${vcfGTSuffix} | bcftools query -f'[%CHROM-%POS-%REF-%ALT\t%SAMPLE\t%GT\t%AD\t%QD\t%AC\t%AN\n]' -i'AC>0 & GT="alt"' | grep -f ${variantsTempFile} | sort -k4,4 > ${outputFile}.gt.tsv
        
        # extract annotations
        cmd="bcftools view -r ${position} ${vcfAnnoPrefix}${fileID}${vcfAnnoSuffix} | bcftools query -f'%CHROM-%POS-%REF-%ALT\t%"${annotationsQuery}"\n' | grep -f ${variantsTempFile} | sort -k4,4 > ${outputFile}.anno.tsv"
        # echo $cmd
        echo $cmd > ${outputFile}.sh
        bash ${outputFile}.sh

        # merge annotations and genotypes
        join -j 1 -a1 <(sort -k1,1 ${outputFile}.gt.tsv) <(sort -k1,1 ${outputFile}.anno.tsv) | sed 's/ /\t/g' >> ${outputFile}
      fi
    else
      if [[ ${addHeader} == "T" ]]; then
        echo -e "variant\t${annotationsHeader}" > ${outputFile}
        addHeader=F
      fi
      # for summary counts without genotypes
      cmd="bcftools query -r ${position} -f '%CHROM-%POS-%REF-%ALT\t%"${annotationsQuery}"\n' ${vcfAnnoPrefix}${fileID}${vcfAnnoSuffix} | grep -f ${variantsTempFile} | sort -k2,2g >> ${outputFile}"
      # echo $cmd
      echo $cmd > ${outputFile}.sh
      bash ${outputFile}.sh
    fi

  done < <(tail -n+2 ${chrGeneList} | grep -v "^#" | grep -E -v "^\s+$" | head -n${topK})
  
  # rm ${variantsTempFile} ${outputFile}.sh 
  # if [[ ${fullGenotype} == "T" && ${vcfGTPrefix} != "" ]]; then
  #   rm ${outputFile}.gt.tsv ${outputFile}.anno.tsv
  # fi
  
}

main "$@"

