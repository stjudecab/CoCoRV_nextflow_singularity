module load nextflow/23.04.1
module load java/openjdk-11
export PYSPARK_SUBMIT_ARGS="--driver-memory 40G pyspark-shell"
inputConfig=/home/stithi/CoCoRV_pipeline_bitbucket/cocorv_pipeline/CoCoRV_nextflow_singularity/example/input.SJLIFE-CNS.GRCh37.txt
nextflow run CoCoRVPipeline.nf -c ${inputConfig} -profile cluster_singularity -with-report ${inputConfig}.report.html -with-timeline ${inputConfig}.timeline.html -with-dag ${inputConfig}.dag.png -resume
