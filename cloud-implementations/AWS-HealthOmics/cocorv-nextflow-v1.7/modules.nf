process coverageIntersect {
  tag "${caseBed}"
  publishDir "/mnt/workflow", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

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
  //tag "${caseVCFPrefix}_${chr}"
  publishDir "/mnt/workflow/vcf_vqsr_normalizedQC", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

  input:
    path case_vcf_file
    path refFASTA

  output:
    path("*.biallelic.leftnorm.ABCheck.vcf.gz")
    path("*.biallelic.leftnorm.ABCheck.vcf.gz.tbi")

  script:
  """
  bash /opt/cocorv/utilities/vcfQCAndNormalize_aws.sh ${case_vcf_file} ${refFASTA}
  """
}


process annotate {
  tag "${chr}"
  publishDir "/mnt/workflow/annotation", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
    path(normalizedQCedVCFFile)
    path(indexFile)
    val build

  output:
    val("${chr}")
    path("${chr}.annotated.vcf.gz")
    path("${chr}.annotated.vcf.gz.tbi")

  script:
  chr = normalizedQCedVCFFile.simpleName
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
  bash /opt/cocorv/utilities/annotate_docker.sh ${normalizedQCedVCFFile} /opt/annovar ${refbuild} \${outputPrefix} ${params.VCFAnno} ${params.toml} ${params.protocol} ${params.operation}

  if [[ ${params.addVEP} == "T" ]]; then
    module load perl/5.28.1
    module load htslib/1.10.2
    module load bcftools/1.15.1
    module load samtools/1.10
    bash /opt/cocorv/utilities/annotateVEPWithOptions.sh ${chr}.annotated.annovar.vcf.gz ${build} ${chr}.annotated ${params.reference} ${params.vepFolder} ${params.cache} ${params.lofteeFolder} ${params.lofteeDataFolder} ${params.caddSNV} ${params.caddIndel} ${params.spliceAISNV} ${params.spliceAIIndel} ${params.perlThread} ${params.AM} ${params.REVEL} 1 ${params.VEPAnnotations}
  fi
  """
}
 
process skipAnnotation {
  tag "${chr}"
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

  input:
    path(normalizedQCedVCFFile)
    path(indexFile)

  output:
    val("${chr}")
    path("${chr}.annotated.vcf.gz")
    path("${chr}.annotated.vcf.gz.tbi")

  script:
  if (chr == "NA") {
    annotated = params.caseAnnotatedVCFFolder + params.caseAnnotatedVCFPrefix + params.caseAnnotatedVCFSuffix
  } else {
    annotated = params.caseAnnotatedVCFFolder + params.caseAnnotatedVCFPrefix + chr + params.caseAnnotatedVCFSuffix
  }
  """
  ln -s ${annotated} ${chr}.annotated.vcf.gz
  """
}


process caseGenotypeGDS {
  tag "${chr}"
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-r:v4'

  cpus params.cpus
  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  publishDir "/mnt/workflow/vcf_vqsr_normalizedQC", mode: 'copy'

  input: 
    path(normalizedQCedVCFFile)
    path(indexFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds")

  script:
  chr = normalizedQCedVCFFile.simpleName
  """
  Rscript /opt/cocorv/utilities/vcf2gds.R ${normalizedQCedVCFFile} ${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds ${params.cpus}
  """
}

process caseAnnotationGDS {
  tag "${chr}"
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-r:v4'

  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  publishDir "/mnt/workflow/annotation", mode: 'copy'

  input: 
    val(chr)
    path(annotatedFile)
    path(indexFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz.gds")

  script:
  """
  Rscript /opt/cocorv/utilities/vcf2gds.R ${annotatedFile} ${chr}.annotated.vcf.gz.gds 1
  """
}

process extractGnomADPositions {
  tag "${chr}"
  publishDir "/mnt/workflow/gnomADPosition", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

  input: 
    path(normalizedQCedVCFFile)
    path(indexFile)
    path(gnomADPCPosition)

  output: 
    path "${chr}.extracted.vcf.gz"

  script:
  chr = normalizedQCedVCFFile.simpleName
  """
  bcftools view -R ${gnomADPCPosition} -Oz -o ${chr}.extracted.vcf.gz ${normalizedQCedVCFFile}
  """
}

process mergeExtractedPositions {
  publishDir "/mnt/workflow/gnomADPosition", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

  input: 
    path extractedVCFFiles

  output: 
    path("all.extracted.vcf.bgz")

  script:
  """
  bcftools concat -Oz -o "all.extracted.vcf.bgz" ${extractedVCFFiles}
  """
}

process RFPrediction {
  publishDir "/mnt/workflow/gnomADPosition", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

  input: 
    path VCFForPrediction
    path loadingPath
    path rfModelPath
    val build
    val threshold

  output: 
    path "PC.population.output.gz"
    path "casePopulation.txt"

  script:
  """
  bash /opt/cocorv/utilities/gnomADPCAndAncestry_docker.sh /opt/cocorv ${loadingPath} ${rfModelPath} ${VCFForPrediction} ${build} "PC.population.output.gz" ${threshold} "casePopulation.txt"
  """
}

process CoCoRV {
  tag "$chr"
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-r:v4'

  memory { 64.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  input: 
    tuple val(chr), path(caseGenotypeGDS), path(caseAnnoGDS) 
    path intersectBed
    path ancestryFile
    path controlAnno_files
    path controlCount_files
    val build
    path ACANConfig
    path caseSample
    val controlAnnoPrefix
    val controlAnnoSuffix
    val controlCountPrefix
    val controlCountSuffix
    path variantExcludeFile
    path highLDVariantFile

  output: 
    path("${chr}.association.tsv") 
    path("${chr}.case.group")
    path("${chr}.control.group")

  script:
  if (chr == "NA") {
    controlAnnotated = controlAnnoPrefix + controlAnnoSuffix
    controlCount = controlCountPrefix + controlCountSuffix
  } else {
    controlAnnotated = controlAnnoPrefix + chr + controlAnnoSuffix
    controlCount = controlCountPrefix + chr + controlCountSuffix
  }
  """
  otherOptions=""
  if [[ "${params.CoCoRVOptions}" != "NA" ]]; then
    otherOptions="${params.CoCoRVOptions}"
  fi
  if [[ "${params.groupFunctionConfig}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantGroupCustom ${params.groupFunctionConfig}"
  fi
  if [[ "${params.extraParamJason}" != "NA" ]]; then
    otherOptions="\${otherOptions} --extraParamJason ${params.extraParamJason}"
  fi
  if [[ "${params.annotationUsed}" != "NA" ]]; then
    otherOptions="\${otherOptions} --annotationUsed ${params.annotationUsed}"
  fi
  if [[ "${params.gnomADVersion}" == "v2exome" ]]; then
    otherOptions="\${otherOptions} --highLDVariantFile ${highLDVariantFile}"
  fi
  if [[ "${params.gnomADVersion}" == "v2exome" || "${params.gnomADVersion}" == "v4exome" ]]; then
    otherOptions="\${otherOptions} --variantExclude ${variantExcludeFile}"
  fi
  if [[ "${params.variantInclude}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantInclude ${params.variantInclude}"
  fi
 
  Rscript /opt/cocorv/utilities/CoCoRV_wrapper.R \
      --sampleList ${caseSample} \
      --outputPrefix ${chr} \
      --AFMax ${params.AFMax} \
      --bed ${intersectBed} \
      --variantMissing ${params.variantMissing} \
      --variantGroup ${params.variantGroup} \
      --removeStar \
      --ACANConfig ${ACANConfig} \
      --caseGroup ${ancestryFile} \
      --minREVEL ${params.REVELThreshold} \
      --checkHighLDInControl \
      --pLDControl ${params.pLDControl} \
      --fullCaseGenotype \
      --reference ${build} \
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
  publishDir "/mnt/workflow/CoCoRV", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-r:v4'

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
  publishDir "/mnt/workflow/CoCoRV", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-r:v4'

  memory { 16.GB * task.attempt }
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
  Rscript /opt/cocorv/utilities/QQPlotAndFDR.R "association.tsv" \
           "association.tsv.dominant.nRep1000" --setID gene \
      --outColumns gene --n 1000 \
      --pattern "case.*Mutation.*_DOM\$|control.*Mutation.*_DOM\$" \
      --FDR
  """
}

process postCheck {
  publishDir "/mnt/workflow/CoCoRV", mode: 'copy'
  container '211125359574.dkr.ecr.us-east-1.amazonaws.com/cocorv-nextflow-python:v6'

  memory { 10.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  input: 
    path associationResult
    val topK
    val caseControl
    val reference
    path caseSample
    path normalizedQCedVCFFiles
    path vcfIndexFiles
    path annotatedFiles
    path annotateIndexFiles
    path caseVariants
    path controlVariants

  output:
    path "*.variants.tsv"

  script:
  """
  bash /opt/cocorv/utilities/checkFPGenes.sh /opt/cocorv/ ${associationResult} ${topK} ${caseControl} ${reference} ${caseSample}
  """
}
