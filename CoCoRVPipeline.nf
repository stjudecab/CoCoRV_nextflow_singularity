// assume tools required are installed
// assume the correct version of scikit-learn is installed

// DSL 2 syntax
nextflow.enable.dsl = 2

// Help Section
log.info """
Usage:
    nextflow run CoCoRVPipeline.nf -c inputConfig -profile profileName
Input:
    * -c: Path of input config file.
    * -profile: Specify which profile to use when executing nextflow pipeline. Available options are: "local", "cluster", "conda_gnomadv4", "cluster_singularity_lsf"
"""

// Import modules
include {coverageIntersect;
  normalizeQC;
  annotate;
  skipAnnotation;
  caseGenotypeGDS;
  caseAnnotationGDS;
  skipGenotypeGDS;
  skipAnnotationGDS;
  extractGnomADPositions;
  mergeExtractedPositions;
  RFPrediction;
  CoCoRV;
  mergeCoCoRVResults;
  QQPlotAndFDR;
  postCheck} from './modules.nf'

// main pipeline
workflow {
  // coverage  
  if (params.caseBed == "NA") {
    intersectChannel = Channel.value(params.controlBed)
  } else {
    coverageIntersect(params.caseBed, params.controlBed)
    intersectChannel = coverageIntersect.out
  }

  // normalize and QC
  chromosomes = params.chrSet.split("\\s+")
  chromChannel = Channel.fromList(Arrays.asList(chromosomes))
  normalizeQC(params.caseVCFPrefix, chromChannel, params.caseVCFSuffix)

  // annotate
  if (params.caseGenotypeGDSPrefix == "NA" && params.caseAnnotationGDSPrefix == "NA") {
    // annotate
    if (params.caseAnnotatedVCFPrefix == "NA") {
      annotate(normalizeQC.out, params.reference)
      annotateChannel = annotate.out
    } else {
      skipAnnotation(normalizeQC.out)
      annotateChannel = skipAnnotation.out
    }

    // case genoypte vcf to gds
    caseGenotypeGDS(normalizeQC.out)
    caseGenotypeGDSChannel = caseGenotypeGDS.out

    // case annotation to gds
    caseAnnotationGDS(annotateChannel)
    caseAnnotationGDSChannel = caseAnnotationGDS.out
  }
  else {
    //skip annotation and GDS conversion
    skipGenotypeGDS(normalizeQC.out)
    caseGenotypeGDSChannel = skipGenotypeGDS.out

    skipAnnotationGDS(skipGenotypeGDS.out)
    caseAnnotationGDSChannel = skipAnnotationGDS.out
  }

  // run gnomAD based population prediction
  if (params.casePopulation == "NA") {
    // extract gnomAD positions
    extractGnomADPositions(normalizeQC.out)

    // merge extracted gnomAD positions
    mergeExtractedPositions(extractGnomADPositions.out[0].collect(),
                            extractGnomADPositions.out[1].collect())

    RFPrediction(mergeExtractedPositions.out)
    populationChannel = RFPrediction.out[1]
  } else {
    populationChannel = Channel.value(params.casePopulation)
  }

  // run CoCoRV
  // RFPrediction.out.view()
  CoCoRV(caseGenotypeGDSChannel.join(caseAnnotationGDSChannel), 
    intersectChannel,
    populationChannel)

  // merge CoCoRV results
  mergeCoCoRVResults(CoCoRV.out[0].collect(), CoCoRV.out[1].collect(), 
    CoCoRV.out[2].collect())

  // QQ plot and FDR
  QQPlotAndFDR(mergeCoCoRVResults.out[0], mergeCoCoRVResults.out[1], mergeCoCoRVResults.out[2])
  postCheck(mergeCoCoRVResults.out[0], params.topK, params.caseControl)
}

