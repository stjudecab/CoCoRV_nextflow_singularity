process coverageIntersect {
  tag "${caseBed}"
  publishDir "${params.outputRoot}", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v1'

  input:
    path caseBed
    path controlBed

  output:
    path "intersect.coverage10x.bed.gz"

  script:
  """
  bedtools intersect -sorted -a ${caseBed} -b ${controlBed} > \
        "intersect.coverage10x.bed.gz"
  """
}

process normalizeQC {
  tag "${caseVCFPrefix}_${chr}"
  publishDir "${params.outputRoot}/vcf_vqsr_normalizedQC", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v1'

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
    vcfFile=caseVCFPrefix + caseVCFSuffix
  } else {
    vcfFile=caseVCFPrefix + chr + caseVCFSuffix
  }
  """
  outputPrefix=${chr}
  ${params.CoCoRVFolder}/utilities/vcfQCAndNormalize.sh ${vcfFile} \${outputPrefix} ${params.refFASTA}
  """
}


process annotate {
  tag "${chr}"
  publishDir "${params.outputRoot}/annotation", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v1'

  memory { 20.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)
    val build

  output:
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz")

  script:
  refbuild = null
  if (build == "GRCh37") {
    refbuild="hg19"  
  }
  else {
    if (build == "GRCh38") {
     refbuild="hg38"   
    }
  }
  """
  if [[ ${params.addVEP} != "T" ]]; then
    outputPrefix="${chr}.annotated"
  else 
    outputPrefix="${chr}.annotated.annovar"
  fi
  bash ${params.CoCoRVFolder}/utilities/annotate.sh ${normalizedQCedVCFFile} ${params.annovarFolder} ${refbuild} \${outputPrefix} ${params.VCFAnno} ${params.toml} ${params.protocol} ${params.operation}

  if [[ ${params.addVEP} == "T" ]]; then
    module load perl/5.28.1
    module load htslib/1.10.2
    module load bcftools/1.15.1
    module load samtools/1.10
    bash ${params.CoCoRVFolder}/utilities/annotateVEPWithOptions.sh ${chr}.annotated.annovar.vcf.gz ${build} ${chr}.annotated ${params.reference} ${params.vepFolder} ${params.cache} ${params.lofteeFolder} ${params.lofteeDataFolder} ${params.caddSNV} ${params.caddIndel} ${params.spliceAISNV} ${params.spliceAIIndel} ${params.perlThread} ${params.AM} ${params.REVEL} 1 ${params.VEPAnnotations}
  fi
  """
}
 
process skipAnnotation {
  tag "${chr}"
  container 'stithi/cocorv-nextflow-python:v1'

  input:
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output:
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz")

  script:
  if (chr == "NA") {
    annotated=params.caseAnnotatedVCFPrefix+params.caseAnnotatedVCFSuffix
  } else {
    annotated=params.caseAnnotatedVCFPrefix+chr+params.caseAnnotatedVCFSuffix
  }
  """
  ln -s ${annotated} ${chr}.annotated.vcf.gz
  """
}


process caseGenotypeGDS {
  tag "${chr}"
  container 'stithi/cocorv-nextflow-r:v1'

  cpus 16
  memory { 20.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  publishDir "${params.outputRoot}/vcf_vqsr_normalizedQC", mode: 'copy'

  input: 
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds")

  script:
  """
  Rscript ${params.CoCoRVFolder}/utilities/vcf2gds.R ${normalizedQCedVCFFile} ${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds 16
  """
}

process caseAnnotationGDS {
  tag "${chr}"
  container 'stithi/cocorv-nextflow-r:v1'

  memory { 20.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  publishDir "${params.outputRoot}/annotation", mode: 'copy'

  input: 
    tuple val(chr), path(annotatedFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz.gds")

  script:
  """
  Rscript ${params.CoCoRVFolder}/utilities/vcf2gds.R ${annotatedFile} ${chr}.annotated.vcf.gz.gds 1
  """
}

process extractGnomADPositions {
  tag "${chr}"
  publishDir "${params.outputRoot}/gnomADPosition", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v1'

  input: 
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output: 
    path "${chr}.extracted.vcf.gz"

  script:
  """
  bcftools view -R ${params.gnomADPCPosition} -Oz -o ${chr}.extracted.vcf.gz ${normalizedQCedVCFFile}
  """
}

process mergeExtractedPositions {
  publishDir "${params.outputRoot}/gnomADPosition", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v1'

  input: 
    path extractedVCFFile

  output: 
    path("all.extracted.vcf.bgz")

  script:
  """
  bcftools concat -Oz -o "all.extracted.vcf.bgz" ${extractedVCFFile}
  """
}

process RFPrediction {
  publishDir "${params.outputRoot}/gnomADPosition", mode: 'copy'
  container 'stithi/cocorv-nextflow-python:v1'

  input: 
    path VCFForPrediction

  output: 
    path "PC.population.output.gz"
    path "casePopulation.txt"

  script:
  """
  bash ${params.CoCoRVFolder}/utilities/gnomADPCAndAncestry.sh ${params.CoCoRVFolder} ${params.loadingPath} ${params.rfModelPath} ${VCFForPrediction} ${params.build} "PC.population.output.gz" ${params.threshold} "casePopulation.txt"
  """
}

process CoCoRV {
  tag "$chr"
  container 'stithi/cocorv-nextflow-r:v1'

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
  if (chr == "NA") {
    controlAnnotated=params.controlAnnoPrefix+params.controlAnnoSuffix
    controlCount=params.controlCountPrefix+params.controlCountSuffix 
  } else {
    controlAnnotated=params.controlAnnoPrefix+chr+params.controlAnnoSuffix
    controlCount=params.controlCountPrefix+chr+params.controlCountSuffix 
  }
  """
  otherOptions="${params.CoCoRVOptions}"
  if [[ "${params.groupFunctionConfig}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantGroupCustom ${params.groupFunctionConfig}"
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
  if [[ "${params.reference}" != "GRCh37" ]]; then
    otherOptions="\${otherOptions} --reference ${params.reference}"
  fi
  if [[ "${params.gnomADVersion}" != "v2exome" ]]; then
    otherOptions="\${otherOptions} --gnomADVersion ${params.gnomADVersion}"
  fi
  if [[ "${params.variantInclude}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantInclude ${params.variantInclude}"
  fi

  Rscript ${params.CoCoRVFolder}/utilities/CoCoRV_wrapper.R \
      --sampleList ${params.caseSample} \
      --outputPrefix ${chr} \
      --AFMax ${params.AFMax} \
      --bed ${intersectBed} \
      --variantMissing ${params.variantMissing} \
      --variantGroup ${params.variantGroup} \
      --removeStar \
      --ACANConfig ${params.ACANConfig} \
      --caseGroup ${ancestryFile} \
      --minREVEL ${params.REVELThreshold} \
      --variantExcludeFile ${params.variantExclude}  \
      --checkHighLDInControl \
      --pLDControl ${params.pLDControl} \
      --fullCaseGenotype \
      --controlAnnoGDSFile ${controlAnnotated} \
      --caseAnnoGDSFile ${caseAnnoGDS} \
      \${otherOptions} \
      ${controlCount} \
      ${caseGenotypeGDS}
  """
}

process mergeCoCoRVResults {
  publishDir "${params.outputRoot}/CoCoRV", mode: 'copy'
  container 'stithi/cocorv-nextflow-r:v1'

  input:
    path associationResult
    path caseVariants
    path controlVariants

  output:
    tuple path("association.tsv"), 
     path("kept.variants.case.txt"), path("kept.variants.control.txt")

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
  cat "association.tsv.tmp" | sort -gk2,2 > "association.tsv"
  """
}

process QQPlotAndFDR {
  publishDir "${params.outputRoot}/CoCoRV", mode: 'copy'
  container 'stithi/cocorv-nextflow-r:v1'

  memory { 10.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  input: 
    tuple path("association.tsv"), 
     path("kept.variants.case.txt"), path("kept.variants.control.txt")

  output:
    path "association.tsv.dominant.nRep1000*"

  script:
  """
  Rscript ${params.CoCoRVFolder}/utilities/QQPlotAndFDR.R "association.tsv" \
           "association.tsv.dominant.nRep1000" --setID gene \
      --outColumns gene --n 1000 \
      --pattern "case.*Mutation.*_DOM\$|control.*Mutation.*_DOM\$" \
      --FDR
  """
}


