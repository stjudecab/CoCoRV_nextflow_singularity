process coverageIntersect {
  tag "${caseBed}"
  publishDir "${params.outputRoot}", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v5'

  input:
    path caseBed
    path controlBed

  output:
    path "intersect.coverage10x.bed.gz"

  script:
  """
  bedtools intersect -a ${caseBed} -b ${controlBed} | gzip > \
        "intersect.coverage10x.bed.gz"
  """
}

process normalizeQC {
  tag "${caseVCFPrefix}_${chr}"
  publishDir "${params.outputRoot}/vcf_vqsr_normalizedQC", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v5'

  input:
    val caseVCFPrefix
    val chr
    val caseVCFSuffix

  output:
    tuple val("${chr}"),
          path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz"),
          path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.tbi")

  script:
  if (chr == "NA") {
    vcfFile = caseVCFPrefix + caseVCFSuffix
  } else {
    vcfFile = caseVCFPrefix + chr + caseVCFSuffix
  }
  """
  outputPrefix=${chr}
  ${params.CoCoRVFolder}/utilities/vcfQCAndNormalize.sh ${vcfFile} \${outputPrefix} ${params.refFASTA}
  """
}


process annotate {
  tag "${chr}"
  publishDir "${params.outputRoot}/annotation", mode: 'copy', pattern: '*.vcf.gz*'
  publishDir "${params.outputRoot}/annotation/logs", mode: 'copy', overwrite: true, pattern: '*.log'
  container 'stithi/cocorv-nextflow-vep:v3'

  cpus params.vepThreads
  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)
    val reference

  output:
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz"),
          path("${chr}.annotated.vcf.gz.tbi")
    path("${chr}.annotated_vep.log"), optional: true

  script:
  refbuild = null
  if (reference == "GRCh37") {
    refbuild="hg19"  
  }
  else {
    if (reference == "GRCh38") {
     refbuild="hg38"   
    }
  }
  """
  if [[ ${params.annotationTool} == "ANNOVAR" ]]; then
    outputPrefix="${chr}.annotated"
    bash ${params.CoCoRVFolder}/utilities/annotate_docker.sh ${normalizedQCedVCFFile} ${params.annovarFolder} ${refbuild} \${outputPrefix} ${params.VCFAnno} ${params.toml} ${params.protocol} ${params.operation}
  fi

  if [[ ${params.annotationTool} == "VEP" ]]; then
    outputPrefix="${chr}.annotated"
    bash -x ${params.CoCoRVFolder}/utilities/annotateVEPWithOptions_docker_no_mane_v3.sh ${normalizedQCedVCFFile} ${reference} ${chr}.annotated ${params.refFASTA} ${params.lofteeFolder} ${params.lofteeDataFolder} ${params.caddSNV} ${params.caddIndel} ${params.spliceAISNV} ${params.spliceAIIndel} ${params.AM} ${params.REVEL} ${params.vepThreads} ${params.VEPAnnotations} ${params.VEPCACHE}
  fi
  
  if [[ ${params.annotationTool} == "ANNOVAR_VEP" ]]; then
    outputPrefix="${chr}.annotated.annovar"
    bash ${params.CoCoRVFolder}/utilities/annotate_docker.sh ${normalizedQCedVCFFile} ${params.annovarFolder} ${refbuild} \${outputPrefix} ${params.VCFAnno} ${params.toml} ${params.protocol} ${params.operation}

    bash ${params.CoCoRVFolder}/utilities/annotateVEPWithOptions_docker_no_mane_v3.sh ${chr}.annotated.annovar.vcf.gz ${reference} ${chr}.annotated ${params.refFASTA} ${params.lofteeFolder} ${params.lofteeDataFolder} ${params.caddSNV} ${params.caddIndel} ${params.spliceAISNV} ${params.spliceAIIndel} ${params.AM} ${params.REVEL} ${params.vepThreads} ${params.VEPAnnotations} ${params.VEPCACHE}

   fi

  """
}
 
process skipAnnotation {
  tag "${chr}"
  publishDir "${params.outputRoot}/annotation", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v5'

  input:
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output:
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz"),
          path("${chr}.annotated.vcf.gz.tbi")

  script:
  if (chr == "NA") {
    annotated = params.caseAnnotatedVCFPrefix + params.caseAnnotatedVCFSuffix
  } else {
    annotated = params.caseAnnotatedVCFPrefix + chr + params.caseAnnotatedVCFSuffix
  }
  """
  ln -s ${annotated} ${chr}.annotated.vcf.gz
  ln -s ${annotated}.tbi ${chr}.annotated.vcf.gz.tbi
  """
}

process caseGenotypeGDS {
  tag "${chr}"
  container 'stithi/cocorv-nextflow-r:v5'

  cpus params.cpus
  memory { 20.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  //publishDir "${params.outputRoot}/vcf_vqsr_normalizedQC", mode: 'copy'

  input: 
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds")

  script:
  """
  Rscript ${params.CoCoRVFolder}/utilities/vcf2gds.R ${normalizedQCedVCFFile} ${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds ${params.cpus}
  """
}

process caseAnnotationGDS {
  tag "${chr}"
  container 'stithi/cocorv-nextflow-r:v5'

  memory { 20.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  //publishDir "${params.outputRoot}/annotation", mode: 'copy'

  input: 
    tuple val(chr), path(annotatedFile), path(indexFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz.gds")

  script:
  """
  Rscript ${params.CoCoRVFolder}/utilities/vcf2gds.R ${annotatedFile} ${chr}.annotated.vcf.gz.gds 1
  """
}

process skipGenotypeGDS {
  tag "${chr}"
  container 'stithi/cocorv-nextflow-python:v5'

  input:
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds")

  script:
  if (chr == "NA") {
    genotypeGDSFile = params.caseGenotypeGDSPrefix + params.caseGenotypeGDSSuffix
  } else {
    genotypeGDSFile = params.caseGenotypeGDSPrefix + chr + params.caseGenotypeGDSSuffix
  }
  """
  ln -s ${genotypeGDSFile} ${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds
  """
}

process skipAnnotationGDS {
  tag "${chr}"
  container 'stithi/cocorv-nextflow-python:v5'

  input:
    tuple val(chr), path(genotypeGDSFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz.gds")

  script:
  if (chr == "NA") {
    annotationGDSFile = params.caseAnnotationGDSPrefix + params.caseAnnotationGDSSuffix
  } else {
    annotationGDSFile = params.caseAnnotationGDSPrefix + chr + params.caseAnnotationGDSSuffix
  }
  """
  ln -s ${annotationGDSFile} ${chr}.annotated.vcf.gz.gds
  """
}

process extractGnomADPositions {
  tag "${chr}"
  publishDir "${params.outputRoot}/gnomADPosition", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v5'

  input: 
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output: 
    path "${chr}.extracted.vcf.gz"
    path "${chr}.extracted.vcf.gz.tbi"

  script:
  """
  bcftools view -R ${params.gnomADPCPosition} -Oz -o ${chr}.extracted.vcf.gz ${normalizedQCedVCFFile}
  bcftools index -t ${chr}.extracted.vcf.gz
  """
}

process mergeExtractedPositions {
  publishDir "${params.outputRoot}/gnomADPosition", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v5'

  input: 
    path extractedVCFFile
    path extractedVCFFileIndex

  output: 
    path("all.extracted.vcf.bgz")

  script:
  """
  bcftools concat -a -Oz -o "all.extracted.vcf.bgz" ${extractedVCFFile}
  """
}

process RFPrediction {
  publishDir "${params.outputRoot}/gnomADPosition", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v5'

  input: 
    path VCFForPrediction

  output: 
    path "PC.population.output.gz"
    path "casePopulation.txt"

  script:
  """
  bash ${params.CoCoRVFolder}/utilities/gnomADPCAndAncestry_docker.sh ${params.CoCoRVFolder} ${params.loadingPath} ${params.rfModelPath} ${VCFForPrediction} ${params.reference} "PC.population.output.gz" ${params.threshold} "casePopulation.txt"
  """
}

process addSexToGroup {
  publishDir "${params.outputRoot}/gnomADPosition", mode: 'copy'
  container 'stithi/cocorv-nextflow-r:v5'

  input: 
    path casePopulation

  output: 
    path "casePopulationBySex.txt"

  script:
  """
  Rscript ${params.CoCoRVFolder}/utilities/stratifiedBySex.R ${casePopulation} ${params.covariate} "casePopulationBySex.txt"
  """
}

process CoCoRV {
  tag "$chr"
  publishDir "${params.outputRoot}/CoCoRV/byChr", mode: 'copy'
  container 'stithi/cocorv-nextflow-r:v5'

  memory { 100.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 0

  input: 
    tuple val(chr), path(caseGenotypeGDS), path(caseAnnoGDS) 
    path intersectBed
    path ancestryFile

  output: 
    path("${chr}.association.tsv") 
    path("${chr}.case.group")
    path("${chr}.control.group")

  script:
  chrOnly = chr
  start = ""
  end = ""
  if (chr == "NA") {
    // NA to use no chr in the controls
    chrOnly = "" 
  } else if (chr.matches(".*_.*")) {
    // this is useful for shad based case data, such as 1_13004384_121976459
    // for chromosome 1 within the region 13004384:121976459, and will match
    // chromosome 1 for the control data
    parts = chr.split("_")
    chrOnly = parts[0]
    start = parts[1]
    end = parts[2]
  } else {
  }
  controlAnnotated=params.controlAnnoPrefix+chrOnly+params.controlAnnoSuffix
  controlCount=params.controlCountPrefix+chrOnly+params.controlCountSuffix
  """
  if [[ "${start}" != "" ]]; then
    # overlap with the shad region
    checkChr=\$(zcat ${intersectBed} | head -1 | cut -f1)
    if [[ \${checkChr} =~ "chr" ]]; then
      chrString="chr"$chrOnly
    else
      chrString=$chrOnly
    fi
    finalBed="intersect.bed.gz"
    printf "\$chrString\t$start\t$end\n" > shad.bed
    bedtools intersect -a ${intersectBed} -b shad.bed | gzip > \${finalBed}
  else
    finalBed=${intersectBed}
  fi

  otherOptions=""
  if [[ "${params.CoCoRVOptions}" != "NA" ]]; then
    otherOptions="${params.CoCoRVOptions}"
  fi
  if [[ "${params.variantGroupCustom}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantGroupCustom ${params.variantGroupCustom}"
  fi
  if [[ "${params.extraParamJason}" != "NA" ]]; then
    otherOptions="\${otherOptions} --extraParamJason ${params.extraParamJason}"
  fi
  if [[ "${params.annotationUsed}" != "NA" ]]; then
    otherOptions="\${otherOptions} --annotationUsed ${params.annotationUsed}"
  fi
  if [[ "${params.highLDVariantFile}" != "NA" ]]; then
    otherOptions="\${otherOptions} --highLDVariantFile ${params.highLDVariantFile}"
  fi
  if [[ "${params.variantExclude}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantExclude ${params.variantExclude}"
  fi
  if [[ "${params.variantInclude}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantInclude ${params.variantInclude}"
  fi

  Rscript ${params.CoCoRVFolder}/utilities/CoCoRV_wrapper.R \
      --sampleList ${params.caseSample} \
      --outputPrefix ${chr} \
      --AFMax ${params.AFMax} \
      --bed \${finalBed} \
      --variantMissing ${params.variantMissing} \
      --groupColumn ${params.groupColumn} \
      --variantGroup ${params.variantGroup} \
      --removeStar \
      --ACANConfig ${params.ACANConfig} \
      --caseGroup ${ancestryFile} \
      --minREVEL ${params.REVELThreshold} \
      --checkHighLDInControl \
      --pLDControl ${params.pLDControl} \
      --fullCaseGenotype \
      --reference ${params.reference} \
      --gnomADVersion ${params.gnomADVersion} \
      --controlAnnoGDSFile ${controlAnnotated} \
      --caseAnnoGDSFile ${caseAnnoGDS} \
      --batchSize ${params.batchSize} \
      --fileID ${chr} \
      \${otherOptions} \
      ${controlCount} \
      ${caseGenotypeGDS}
  """
}

process mergeCoCoRVResults {
  publishDir "${params.outputRoot}/CoCoRV", mode: 'copy'
  container 'stithi/cocorv-nextflow-r:v5'

  input:
    path associationResult
    path caseVariants
    path controlVariants

  output:
    path("association.tsv")
    path("kept.variants.case.txt")
    path("kept.variants.control.txt")

  script:
  """
  body() {
    IFS= read -r header
    printf '%s\n' "\$header"
    "\$@"
  }

  i=0
  cat ${caseVariants} > "kept.variants.case.txt"
  cat ${controlVariants} > "kept.variants.control.txt"
  for file in ${associationResult}; do
    echo \${file}
    if [[ \${i} == 0 ]]; then
      cat \${file} > "association.tsv.tmp"
      i=1
    else
      tail -n+2 \${file} >> "association.tsv.tmp"
    fi
  done
  cat "association.tsv.tmp" | sort -gk3,3 > "association.tsv"
  """
}

process QQPlotAndFDR {
  publishDir "${params.outputRoot}/CoCoRV", mode: 'copy'
  container 'stithi/cocorv-nextflow-r:v5'

  memory { 10.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  input: 
    path("association.tsv")
    path("kept.variants.case.txt")
    path("kept.variants.control.txt")

  output:
    path "association.tsv.dominant.nRep1000*"

  script:
  """
  Rscript ${params.CoCoRVFolder}/utilities/QQPlotAndFDR.R "association.tsv" \
           "association.tsv.dominant.nRep1000" --setID gene \
      --outColumns gene,P_DOM,OR_DOM,CI_Lower_DOM,CI_Upper_DOM --n 1000 \
      --pattern "case.*Mutation.*_DOM\$|control.*Mutation.*_DOM\$" \
      --FDR
  """
}

process postCheck {
  publishDir "${params.outputRoot}/CoCoRV", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v5'

  memory { 10.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  input: 
    path associationResult
    val topK
    val caseControl

  output:
    path "*.variants.tsv"

  script:
  """
  if [[ ${params.annotationTool} == "ANNOVAR" ]]; then
    bash ${params.CoCoRVFolder}/utilities/checkFPGenes.sh ${params.CoCoRVFolder} ${associationResult} ${topK} ${caseControl} ${params.outputRoot} ${params.reference} ${params.caseSample} "F"
  elif [[ ${params.annotationTool} == "VEP" || ${params.annotationTool} == "ANNOVAR_VEP" ]]; then
    bash ${params.CoCoRVFolder}/utilities/checkFPGenes.sh ${params.CoCoRVFolder} ${associationResult} ${topK} ${caseControl} ${params.outputRoot} ${params.reference} ${params.caseSample} "T"
  fi
  """
}
