module load nextflow/23.04.1
module load java/openjdk-11
export PYSPARK_SUBMIT_ARGS="--driver-memory 40G pyspark-shell"
inputConfig=/home/stithi/CoCoRV_pipeline_bitbucket/cocorv_pipeline/CoCoRV_nextflow_singularity/example/input.CReATe-ALS.gnomAD.v4.1exomes-sex-stratified.txt
nextflow run CoCoRVPipeline.nf -c ${inputConfig} -profile cluster_singularity_lsf -w /scratch_space/stithi/Nextflow_work/v4.1exome -with-report ${inputConfig}.report.html -with-timeline ${inputConfig}.timeline.html -with-dag ${inputConfig}.dag.png -resume
