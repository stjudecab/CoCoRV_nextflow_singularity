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
    * -profile: Specify which profile to use when executing nextflow pipeline. Available options are: "local", "cluster", "conda_gnomADv4"
"""

// Import modules
include {coverageIntersect;
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
  postCheck} from './modules.nf'

// main pipeline
workflow {
  //set the params
  if (params.gnomADVersion == "v4exome") {
    build = "GRCh38"
    refFASTA = "s3://cocorv-resource-files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
    //annovarFolder = "s3://cocorv-annovar-data/"
    
    controlBed = "s3://cocorv-resource-files/gnomADv4.1exome/coverage10x.bed.gz"
    controlAnnoFolder = "s3://cocorv-resource-files/gnomADv4.1exome/annotation/"
    controlAnnoPrefix = "chr"
    controlAnnoSuffix = ".annovar.vep.vcf.gz.gds"
    controlCountFolder = "s3://cocorv-resource-files/gnomADv4.1exome/genotypeCount/"
    controlCountPrefix = "gnomad.exomes.v4.1.sites.chr"
    controlCountSuffix = ".vcf.bgz.gds"
    
    gnomADPCPosition = "s3://cocorv-resource-files/ancestryPrediction/hail_positions.GRCh38.v4.chr.pos.tsv"
    loadingPath = "s3://cocorv-resource-files/ancestryPrediction/gnomad.v4.0.pca_loadings.ht/"
    rfModelPath = "s3://cocorv-resource-files/ancestryPrediction/gnomad.v4.0.RF_fit.onnx"
    threshold = "0.75"

    variantExcludeFile = "s3://cocorv-resource-files/settings-files/gnomAD41WGSExtraExcludeInCodingExcludeTAS2R46.txt.gz"    
  }
  else if (params.gnomADVersion == "v4genome") {
    build = "GRCh38"
    refFASTA = "s3://cocorv-resource-files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
    //annovarFolder = "s3://cocorv-annovar-data/"
    
    controlBed = "s3://cocorv-resource-files/gnomADv4.1genome/coverage10x.bed.gz"
    controlAnnoFolder = "s3://cocorv-resource-files/gnomADv4.1genome/annotation/"
    controlAnnoPrefix = "chr"
    controlAnnoSuffix = ".annovar.vep.vcf.gz.gds"
    controlCountFolder = "s3://cocorv-resource-files/gnomADv4.1genome/genotypeCount/"
    controlCountPrefix = "gnomad.genomes.v4.1.sites.chr"
    controlCountSuffix = ".vcf.bgz.gds"
    
    gnomADPCPosition = "s3://cocorv-resource-files/ancestryPrediction/hail_positions.GRCh38.v4.chr.pos.tsv"
    loadingPath = "s3://cocorv-resource-files/ancestryPrediction/gnomad.v4.0.pca_loadings.ht/"
    rfModelPath = "s3://cocorv-resource-files/ancestryPrediction/gnomad.v4.0.RF_fit.onnx"
    threshold = "0.75"

    variantExcludeFile = "s3://cocorv-resource-files/settings-files/gnomAD41WGSExtraExcludeInCodingExcludeTAS2R46.txt.gz"
  }
  else if (params.gnomADVersion == "v2exome") {
    build = "GRCh37"
    refFASTA = "s3://cocorv-resource-files/GRCh37-lite.fa"
    //annovarFolder = "s3://cocorv-annovar-data/"
    
    controlBed = "s3://cocorv-resource-files/gnomADv2exome/coverage10x.bed.gz"
    controlAnnoFolder = "s3://cocorv-resource-files/gnomADv2exome/annotation/"
    controlAnnoPrefix = "chr"
    controlAnnoSuffix = ".vep.fixed.vcf.gz.gds"
    controlCountFolder = "s3://cocorv-resource-files/gnomADv2exome/genotypeCount/"
    controlCountPrefix = "gnomad.exomes.r2.1.1.sites."
    controlCountSuffix = ".vcf.bgz.gds"
    
    gnomADPCPosition = "s3://cocorv-resource-files/ancestryPrediction/hail_positions.GRCh37.chr.pos.tsv"
    loadingPath = "s3://cocorv-resource-files/ancestryPrediction/gnomad.r2.1.pca_loadings.ht/"
    rfModelPath = "s3://cocorv-resource-files/ancestryPrediction/gnomad.r2.1.RF_fit.onnx"
    threshold = "0.9"

    variantExcludeFile = "s3://cocorv-resource-files/settings-files/gnomAD.exclude.allow.segdup.lcr.v3.txt.gz"  
  }
  highLDVariantFile = "s3://cocorv-resource-files/settings-files/full_vs_gnomAD.p0.05.OR1.ignoreEthnicityInLD.rds"

  // coverage  
  if (params.caseBed == "NA") {
    intersectChannel = Channel.value(controlBed)   
  }
  else {
    coverageIntersect(params.caseBed, controlBed)
    intersectChannel = coverageIntersect.out
  }

  // normalize and QC
  //chromosomes = params.chrSet.split("\\s+")
  //chromChannel = Channel.fromList(Arrays.asList(chromosomes))
  //case_vcf_ch = Channel
  //    .fromPath(params.caseVCFFolder + params.caseVCFPrefix + '*' + params.caseVCFSuffix)
  case_vcf_ch = Channel            
			            .fromPath(params.caseVCFFileList)               
                  .splitText()
                  .map{it.replaceFirst(/\n/,"")}
                  .map{ file(it) }   //map the file path string into file object, then can extract the file information. 
                  
  normalizeQC(case_vcf_ch, refFASTA)

  // annotate
  if (params.caseAnnotatedVCFPrefix == "NA") {
    annotate(normalizeQC.out[0], normalizeQC.out[1], build)
    annotateChannel = annotate.out
  } else {
    skipAnnotation(normalizeQC.out)
    annotateChannel = skipAnnotation.out
  }

  // case genoypte vcf to gds
  caseGenotypeGDS(normalizeQC.out[0], normalizeQC.out[1])

  // case annotation to gds
  caseAnnotationGDS(annotateChannel)

  // run gnomAD based population prediction
  if (params.casePopulation == "NA") {
    // extract gnomAD positions
    extractGnomADPositions(normalizeQC.out[0], normalizeQC.out[1], gnomADPCPosition)

    // merge extracted gnomAD positions
    mergeExtractedPositions(extractGnomADPositions.out.collect())

    RFPrediction(mergeExtractedPositions.out, loadingPath, rfModelPath, build, threshold)
    populationChannel = RFPrediction.out[1]
  } else {
    populationChannel = Channel.value(params.casePopulation)
  }

  // run CoCoRV
  // RFPrediction.out.view()
  controlAnno_ch = Channel
      .fromPath(controlAnnoFolder + controlAnnoPrefix + "*" + controlAnnoSuffix)
  controlCount_ch = Channel
      .fromPath(controlCountFolder + controlCountPrefix + "*" + controlCountSuffix)
  
  CoCoRV(caseGenotypeGDS.out.join(caseAnnotationGDS.out), 
    intersectChannel,
    populationChannel,
    controlAnno_ch.collect(),
    controlCount_ch.collect(),
    build,
    params.ACANConfig,
    params.caseSample,
    controlAnnoPrefix,
    controlAnnoSuffix,
    controlCountPrefix,
    controlCountSuffix,
    variantExcludeFile,
    highLDVariantFile)

  // merge CoCoRV results
  mergeCoCoRVResults(CoCoRV.out[0].collect(), CoCoRV.out[1].collect(), 
    CoCoRV.out[2].collect())

  // QQ plot and FDR
  QQPlotAndFDR(mergeCoCoRVResults.out[0], mergeCoCoRVResults.out[1], mergeCoCoRVResults.out[2])
  postCheck(mergeCoCoRVResults.out[0], params.topK, params.caseControl, build, params.caseSample, 
    normalizeQC.out[1].collect(), normalizeQC.out[2].collect(), annotate.out[1].collect(), annotate.out[2].collect(), 
    CoCoRV.out[1].collect(), CoCoRV.out[2].collect())
}

