process splitJointVCF {
    tag "${caseJointVCF}_${chr}"
    label 'process_single'
    container 'stithi/cocorv-nextflow-python:v7'

    input:
    path caseJointVCF
    path caseJointVCFtbi
    val chr

    output:
    tuple val("${chr}"), path("chr${chr}.vcf.gz")

    script:
    """
    if [[ ${params.build} == "GRCh37" ]]; then
        bcftools view -r ${chr} -Oz -o chr${chr}.vcf.gz ${caseJointVCF}
    else
        bcftools view -r chr${chr} -Oz -o chr${chr}.vcf.gz ${caseJointVCF}
    fi
    """
}

process coverageIntersect {
    tag "${caseBed}"
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy'
    container 'stithi/cocorv-nextflow-python:v7'

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
    tag "${chr}"
    label 'process_single'
    publishDir "${params.outdir}/vcf_vqsr_normalizedQC", mode: 'copy'
    container 'stithi/cocorv-nextflow-python:v7'

    input:
    tuple val(chr), path(vcfFile)
    path reference
    path referenceFai
    path referenceGzi

    output:
    val("${chr}"), emit: chr
    path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz"), emit: normalizedQCedVCFFile
    path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.tbi"), emit: normalizedQCedVCFFileIndex

    script:
    """
    outputPrefix=${chr}
    ${params.CoCoRVFolder}/utilities/vcfQCAndNormalize.sh ${vcfFile} \${outputPrefix} ${reference}
    """
}


process annotate {
    tag "${chr}"
    label 'process_medium'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    container 'stithi/cocorv-nextflow-vep:v3'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries 5

    input:
    val(chr)
    path(normalizedQCedVCFFile)
    path(indexFile)
    val build
    path annovarFolder
    path vepFolder
    path reference

    output:
    val("${chr}"), emit: chr
    path("${chr}.annotated.vcf.gz"), emit: annotatedFile
    path("${chr}.annotated.vcf.gz.tbi"), emit: annotatedFileIndex

    script:
    refbuild = null
    caddSNV = null
    caddIndel = null
    spliceAISNV = null
    spliceAIIndel = null
    AM = null
    REVEL = null
    VEPCACHE = null
    vepThreads = 6

    if (build == "GRCh37") {
        refbuild="hg19"
        lofteeFolder = vepFolder + "/other_data/loftee/loftee"
        lofteeDataFolder = vepFolder + "/other_data/loftee/data"
        caddSNV = vepFolder + "/other_data/CADD/hg19/v1.7/whole_genome_SNVs.tsv.gz"
        caddIndel = vepFolder + "/other_data/CADD/hg19/v1.7/gnomad.genomes-exomes.r4.0.indel.tsv.gz"
        spliceAISNV = vepFolder + "/other_data/SpliceAI/spliceai_scores.raw.snv.hg19.vcf.gz"
        spliceAIIndel = vepFolder + "/other_data/SpliceAI/spliceai_scores.raw.indel.hg19.vcf.gz"
        AM = vepFolder + "/other_data/AlphaMissense"
        REVEL = vepFolder + "/other_data/REVEL"
        VEPCACHE = vepFolder + "/ensembl-vep/cache"
    }
    else if (build == "GRCh38") {
        refbuild="hg38"
        lofteeFolder = vepFolder + "/other_data/loftee/loftee"
        lofteeDataFolder = vepFolder + "/other_data/loftee/data"
        caddSNV = vepFolder + "/other_data/CADD/hg38/v1.7/whole_genome_SNVs.tsv.gz"
        caddIndel = vepFolder + "/other_data/CADD/hg38/v1.7/gnomad.genomes-exomes.r4.0.indel.tsv.gz"
        spliceAISNV = vepFolder + "/other_data/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
        spliceAIIndel = vepFolder + "/other_data/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
        AM = vepFolder + "/other_data/AlphaMissense"
        REVEL = vepFolder + "/other_data/REVEL"
        VEPCACHE = vepFolder + "/ensembl-vep/cache_hg38"
    }

    """
    if [[ ${params.addVEP} != "T" ]]; then
        outputPrefix="${chr}.annotated"
    else
        outputPrefix="${chr}.annotated.annovar"
    fi
    bash ${params.CoCoRVFolder}/utilities/annotate_docker.sh ${normalizedQCedVCFFile} ${annovarFolder} ${refbuild} \${outputPrefix} ${params.VCFAnno} ${params.toml} ${params.protocol} ${params.operation}

    if [[ ${params.addVEP} == "T" ]]; then
        bash ${params.CoCoRVFolder}/utilities/annotateVEPWithOptions_docker_no_mane_v3.sh ${chr}.annotated.annovar.vcf.gz ${build} ${chr}.annotated ${reference} ${lofteeFolder} ${lofteeDataFolder} ${caddSNV} ${caddIndel} ${spliceAISNV} ${spliceAIIndel} ${AM} ${REVEL} ${vepThreads} ${params.VEPAnnotations} ${VEPCACHE}
    fi
    """
}

process skipAnnotation {
    tag "${chr}"
    label 'process_single'
    publishDir "${params.outdir}/annotation", mode: 'copy'
    container 'stithi/cocorv-nextflow-python:v7'

    input:
    tuple val(chr), path(annotated), path(annotatedTbi)

    output:
    val("${chr}"), emit: chr
    path("${chr}.annotated.vcf.gz"), emit: annotatedFile
    path("${chr}.annotated.vcf.gz.tbi"), emit: annotatedFileIndex

    script:
    """
    """
}

process caseGenotypeGDS {
    tag "${chr}"
    label 'process_medium'
    container 'stithi/cocorv-nextflow-r:v5'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries 5

    publishDir "${params.outdir}/vcf_vqsr_normalizedQC", mode: 'copy'

    input:
    val(chr)
    path(normalizedQCedVCFFile)
    path(indexFile)

    output:
    tuple val("${chr}"),
    path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds")

    script:
    """
    Rscript ${params.CoCoRVFolder}/utilities/vcf2gds.R ${normalizedQCedVCFFile} ${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds 4
    """
}

process caseAnnotationGDS {
    tag "${chr}"
    label 'process_medium'
    container 'stithi/cocorv-nextflow-r:v5'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries 1

    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    val(chr)
    path(annotatedFile)
    path(indexFile)

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
    label 'process_single'
    publishDir "${params.outdir}/gnomADPosition", mode: 'copy'
    container 'stithi/cocorv-nextflow-python:v7'

    input:
    val(chr)
    path(normalizedQCedVCFFile)
    path(indexFile)
    path(gnomADPCPosition)

    output:
    path "${chr}.extracted.vcf.gz"

    script:
    """
    bcftools view -R ${gnomADPCPosition} -Oz -o ${chr}.extracted.vcf.gz ${normalizedQCedVCFFile}
    """
}

process mergeExtractedPositions {
    label 'process_single'
    publishDir "${params.outdir}/gnomADPosition", mode: 'copy'
    container 'stithi/cocorv-nextflow-python:v7'

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
    label 'process_single'
    publishDir "${params.outdir}/gnomADPosition", mode: 'copy'
    container 'stithi/cocorv-nextflow-python:v7'

    input:
    path VCFForPrediction
    path loadingPath
    path rfModelPath

    output:
    path "PC.population.output.gz"
    path "casePopulation.txt"

    script:
    threshold = "0.75"
    if (params.build == "GRCh37") {
    threshold = "0.9"
    }
    else if (params.build == "GRCh38") {
    threshold = "0.75"
    }
    """
    bash ${params.CoCoRVFolder}/utilities/gnomADPCAndAncestry_docker.sh ${params.CoCoRVFolder} ${loadingPath} ${rfModelPath} ${VCFForPrediction} ${params.build} "PC.population.output.gz" ${threshold} "casePopulation.txt"
    """
}

process CoCoRV {
    tag "$chr"
    label 'process_high_memory'
    publishDir "${params.outdir}/CoCoRV/byChr", mode: 'copy'
    container 'stithi/cocorv-nextflow-r:v5'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries 0

    input:
    tuple val(chr), path(caseGenotypeGDS), path(caseAnnoGDS), path(controlCount), path(controlAnnotated)
    path intersectBed
    path ancestryFile
    path ACANConfig
    path variantExclude
    path highLDVariantFile
    path caseSample

    output:
    path("${chr}.association.tsv"), emit: association_perChr
    path("${chr}.case.group"), emit: caseVariants_perChr
    path("${chr}.control.group"), emit: controlVariants_perChr

    script:

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
        otherOptions="\${otherOptions} --variantExclude ${variantExclude}"
    fi
    if [[ "${params.variantInclude}" != "NA" ]]; then
        otherOptions="\${otherOptions} --variantInclude ${params.variantInclude}"
    fi

    Rscript ${params.CoCoRVFolder}/utilities/CoCoRV_wrapper.R \
    --sampleList ${caseSample} \
    --outputPrefix ${chr} \
    --AFMax ${params.AFMax} \
    --bed ${intersectBed} \
    --variantMissing ${params.variantMissing} \
    --groupColumn ${params.groupColumn} \
    --variantGroup ${params.variantGroup} \
    --removeStar \
    --ACANConfig ${ACANConfig} \
    --caseGroup ${ancestryFile} \
    --minREVEL ${params.REVELThreshold} \
    --checkHighLDInControl \
    --pLDControl ${params.pLDControl} \
    --fullCaseGenotype \
    --reference ${params.build} \
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
    label 'process_medium'
    publishDir "${params.outdir}/CoCoRV", mode: 'copy'
    container 'stithi/cocorv-nextflow-r:v5'

    input:
    path associationResult
    path caseVariants
    path controlVariants

    output:
    path("association.tsv"), emit: association_res
    path("kept.variants.case.txt"), emit: caseVariants_res
    path("kept.variants.control.txt"), emit: controlVariants_res

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
    label 'process_medium'
    publishDir "${params.outdir}/CoCoRV", mode: 'copy'
    container 'stithi/cocorv-nextflow-r:v5'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path("association.tsv")
    path("kept.variants.case.txt")
    path("kept.variants.control.txt")

    output:
    path "association.tsv.dominant.nRep1000.pdf", emit: qqplot
    path "association.tsv.dominant.nRep1000.fdr.tsv", emit: fdr_res

    script:
    """
    Rscript ${params.CoCoRVFolder}/utilities/QQPlotAndFDR.R "association.tsv" \
        "association.tsv.dominant.nRep1000" --setID gene \
        --outColumns gene --n 1000 \
        --pattern "case.*Mutation.*_DOM\$|control.*Mutation.*_DOM\$" \
        --FDR
    """
}

process postCheck {
    label 'process_high_memory'
    publishDir "${params.outdir}/CoCoRV", mode: 'copy'
    container 'stithi/cocorv-nextflow-python:v7'

    errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path associationResult
    val topK
    val caseControl
    val reference
    path caseSample
    path normalizedQCedVCFFiles
    path normalizedindexFiles
    path annotatedFiles
    path annotateindexFiles
    path caseVariants
    path controlVariants

    output:
    path "*.variants.tsv"

    script:
    """
    bash ${params.CoCoRVFolder}/utilities/checkFPGenes.sh ${params.CoCoRVFolder} ${associationResult} ${topK} ${caseControl} ${params.build} ${caseSample} ${params.addVEP}
    """
}
