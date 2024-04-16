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
    * -profile: Specify which profle to use when executing nextflow pipeline. Availble options are: "local", "cluster", "conda_GRCh37", "conda_GRCh38"
"""

// Import modules
include {coverageIntersect;
  normalizeQC;
  annotateANNOVAR;
  skipANNOVAR;
  caseGenotypeGDS;
  caseAnnotationGDS;
  extractGnomADPositions;
  mergeExtractedPositions;
  RFPrediction;
  CoCoRV;
  mergeCoCoRVResults;
  QQPlotAndFDR} from './modules.nf'

// main pipeline
workflow {
  // coverage  
  if (params.caseBed == "NA") {
    intersectChannel = Channel.fromPath(params.controlBed)   
  }
  else {
    coverageIntersect(params.caseBed, params.controlBed)
    intersectChannel = coverageIntersect.out
  }

  // normalize and QC
  chromosomes = params.chrSet.split("\\s+")
  chromChannel = Channel.fromList(Arrays.asList(chromosomes))
  normalizeQC(params.caseVCFPrefix, chromChannel, params.caseVCFSuffix)

  // annotate
  if (params.caseAnnotatedVCFPrefix == "NA") {
    annotateANNOVAR(normalizeQC.out, params.build)
    annotateChannel = annotateANNOVAR.out
  } else {
    skipANNOVAR(normalizeQC.out)
    annotateChannel = skipANNOVAR.out
  }

  // case genoypte vcf to gds
  caseGenotypeGDS(normalizeQC.out)

  // case annotation to gds
  caseAnnotationGDS(annotateChannel)


  // run gnomAD based population prediction
  if (params.casePopulation == "NA") {
    // extract gnomAD positions
    extractGnomADPositions(normalizeQC.out)

    // merge extracted gnomAD positions
    mergeExtractedPositions(extractGnomADPositions.out.collect())

    RFPrediction(mergeExtractedPositions.out)
    populationChannel = RFPrediction.out[1]
  } else {
    populationChannel = Channel.value(params.casePopulation)
  }

  // run CoCoRV
  // RFPrediction.out.view()
  CoCoRV(caseGenotypeGDS.out.join(caseAnnotationGDS.out), 
    intersectChannel,
    populationChannel)

  // merge CoCoRV results
  mergeCoCoRVResults(CoCoRV.out[0].collect(), CoCoRV.out[1].collect(), 
    CoCoRV.out[2].collect())

  // QQ plot and FDR
  QQPlotAndFDR(mergeCoCoRVResults.out)
}

