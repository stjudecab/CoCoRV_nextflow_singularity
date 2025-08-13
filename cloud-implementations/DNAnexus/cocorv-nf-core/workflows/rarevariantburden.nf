/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rarevariantburden_pipeline'

include {
    splitJointVCF;
    coverageIntersect;
    normalizeQC;
    annotate;
    skipAnnotation;
    caseGenotypeGDS;
    caseAnnotationGDS;
    extractGnomADPositions;
    mergeExtractedPositions;
    RFPrediction;
    CoCoRV;
    mergeCoCoRVResults;
    QQPlotAndFDR;
    postCheck
} from '../modules/local/modules.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RAREVARIANTBURDEN {

    take:
    caseJointVCF // caseJointVCF read in from --caseJointVCF
    caseSample // caseSample read in from --caseSample

    main:

    // coverage
    if (params.caseBed == "NA") {
        intersectChannel = Channel.value(params.controlBed)
    } else {
        coverageIntersect(params.caseBed, params.controlBed)
        intersectChannel = coverageIntersect.out
    }

    // create chromosome channel
    chromosomes = params.chrSet.split("\\s+")
    chromChannel = Channel.fromList(Arrays.asList(chromosomes))

    // split joint VCF by chromosome
    // normalize and QC
    if (params.caseVCFFileList == "NA") {
        caseJointVCFtbi = params.caseJointVCF + ".tbi"
        splitJointVCF(params.caseJointVCF, caseJointVCFtbi, chromChannel)
        normalizeQC(splitJointVCF.out, params.refFASTA, params.refFASTA + ".fai", params.refFASTA + ".gzi")
        normalizeQCChannel = normalizeQC.out
    } else {
        caseVCF_ch = Channel
			            .fromPath(params.caseVCFFileList)
                        .splitCsv(header: true)
                        .map { row -> tuple(row.chr, file(row.vcf)) } // Create a tuple of chr and case VCF file path

        normalizeQC(caseVCF_ch, params.refFASTA, params.refFASTA + ".fai", params.refFASTA + ".gzi")
        normalizeQCChannel = normalizeQC.out
    }

    // annotate
    if (params.caseAnnotatedVCFFileList == "NA") {
        
        annotate(normalizeQCChannel, params.build, params.annovarFolder, params.vepFolder, params.refFASTA)
        annotateChannel = annotate.out
    } else {
        annotate_ch = Channel
			                    .fromPath(params.caseAnnotatedVCFFileList)
                                .splitCsv(header: true)
                                .map { row -> tuple(row.chr, file(row.vcf), file(row.index)) } // Create a tuple of chr, annotated VCF file path, index file path

        skipAnnotation(annotate_ch)
        annotateChannel = skipAnnotation.out
    }

    if (params.caseGenotypeGDSFileList == "NA" && params.caseAnnotationGDSFileList == "NA") {
        // case genoypte vcf to gds
        caseGenotypeGDS(normalizeQCChannel)
        caseGenotypeGDSChannel = caseGenotypeGDS.out

        // case annotation to gds
        caseAnnotationGDS(annotateChannel)
        caseAnnotationGDSChannel = caseAnnotationGDS.out
    }
    else {
        //skip annotation and GDS conversion
        caseGenotypeGDSChannel = Channel
			            .fromPath(params.caseGenotypeGDSFileList)
                        .splitCsv(header: true)
                        .map { row -> tuple(row.chr, file(row.gds)) } // Create a tuple of chr and GDS file path


        caseAnnotationGDSChannel = Channel
			            .fromPath(params.caseAnnotationGDSFileList)
                        .splitCsv(header: true)
                        .map { row -> tuple(row.chr, file(row.gds)) } // Create a tuple of chr and GDS file path

    }

    // run gnomAD based population prediction
    if (params.casePopulation == "NA") {
        // extract gnomAD positions
        extractGnomADPositions(normalizeQCChannel, params.gnomADPCPosition)

        // merge extracted gnomAD positions
        mergeExtractedPositions(extractGnomADPositions.out.collect())

        RFPrediction(mergeExtractedPositions.out, params.loadingPath, params.rfModelPath)
        populationChannel = RFPrediction.out[1]
    } else {
        populationChannel = Channel.value(params.casePopulation)
    }

    // run CoCoRV
    // RFPrediction.out.view()

    controlGenotypeGDSChannel = Channel
			            .fromPath(params.controlGenotypeGDSFileList)
                        .splitCsv(header: true)
                        .map { row -> tuple(row.chr, file(row.gds)) } // Create a tuple of chr and GDS file path

    controlAnnotationGDSChannel = Channel
			            .fromPath(params.controlAnnotationGDSFileList)
                        .splitCsv(header: true)
                        .map { row -> tuple(row.chr, file(row.gds)) } // Create a tuple of chr and GDS file path

    controlChannel = controlGenotypeGDSChannel.join(controlAnnotationGDSChannel)
    caseChannel = caseGenotypeGDSChannel.join(caseAnnotationGDSChannel)

    if (params.build == "GRCh37") {
        CoCoRV(caseChannel.join(controlChannel),
            intersectChannel,
            populationChannel,
            params.controlDataFolder + "/stratified_config_gnomad.txt",
            params.controlDataFolder + "/gnomAD.exclude.allow.segdup.lcr.v3.txt.gz",
            params.controlDataFolder + "/full_vs_gnomAD.p0.05.OR1.ignoreEthnicityInLD.rds",
            params.caseSample)
    }
    else if (params.build == "GRCh38") {
        CoCoRV(caseChannel.join(controlChannel),
            intersectChannel,
            populationChannel,
            params.controlDataFolder + "/stratified_config_gnomadV4.asj.txt",
            params.controlDataFolder + "/gnomAD41WGSExtraExcludeInCodingExcludeTAS2R46.txt.gz",
            params.controlDataFolder + "/full_vs_gnomAD.p0.05.OR1.ignoreEthnicityInLD.rds",
            params.caseSample)
    }

    // merge CoCoRV results
    mergeCoCoRVResults(CoCoRV.out.association_perChr.collect(), CoCoRV.out.caseVariants_perChr.collect(),
        CoCoRV.out.controlVariants_perChr.collect())

    // QQ plot and FDR
    QQPlotAndFDR(mergeCoCoRVResults.out.association_res, mergeCoCoRVResults.out.caseVariants_res, mergeCoCoRVResults.out.controlVariants_res)

    //postCheck(mergeCoCoRVResults.out[0], params.topK, params.caseControl)

    postCheck(mergeCoCoRVResults.out.association_res, params.topK, params.caseControl, params.build, params.caseSample,
        normalizeQCChannel.normalizedQCedVCFFile.collect(),
        normalizeQCChannel.normalizedQCedVCFFileIndex.collect(),
        annotateChannel.annotatedFile.collect(),
        annotateChannel.annotatedFileIndex.collect(),
        CoCoRV.out.caseVariants_perChr.collect(), CoCoRV.out.controlVariants_perChr.collect())

    emit:association_res = mergeCoCoRVResults.out.association_res // channel: /path/to/association.tsv
    qqplot               = QQPlotAndFDR.out.qqplot                // channel: /path/to/association.tsv.dominant.nRep1000.pdf

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
