# A nextflow based pipeline for CoCoRV based analysis #

### What is this repository for? ###
This repository includes the nextflow based pipeline for applying CoCoRV to sequencing data, e.g., whole exome sequencing (WES), or whole genome sequencing (WGS), using gnomAD data as control. It supports gnomAD v2 GRCh37 exome data, and gnomAD v4 exome data, and gnomAD v4 genome data as control.    

### Installation ###

#### Nextflow ####
Load nextflow using the following commands:
```bash
module load nextflow/23.04.1
module load java/openjdk-11
```
If you want to install nextflow, please check https://www.nextflow.io/.

#### Singularity ####
Singularity module is already loaded on our HPC cluster. We don't need to load it.

#### Run CoCoRV pipeline ####
Download the git repository using the following command:
```bash
git clone https://github.com/stjudecab/CoCoRV_nextflow_singularity.git
```
Go to the project folder and run main nextflow script *CoCoRVPipeline.nf*
```bash
nextflow run CoCoRVPipeline.nf
```
It will print the help message.

### Examples ###
The CoCoRV nextflow pipeline script, nextflow configuration file, and example run configuration files are provided in this repository.

* *CoCoRVPipeline.nf*: main nextflow script for CoCoRV pipeline.
* *nextflow.config*: contains default configuration for running CoCoRV pipeline.
* *example folder*: contains run specific configuration files. You need to create a configuration file like "input.1KG.GRCh38.gnomAD.v4exomes.txt" or "input.1KG.GRCh38.gnomAD.v4genomes.txt" given in here and update it with case VCF file path and output folder path for your run.

#### GRCh37 using gnmoAD v2 exome data ####
If you have GRCh37 data, you need to use this version. To run the CoCoRV pipeline for GRCh37 using gnomAD v2 exome data, an example run script "testGRCh37.gnomAD.v2exome_1KG_singularity.sh" is given. This test script uses the input configuration file given in here: "example/input.1KG.GRCh37.txt". Here we used 25 test samples from 1000 Genomes Project (build GRCh37) as case data.

To run "testGRCh37.gnomAD.v2exome_1KG_singularity.sh" script, update the "inputConfig" parameter to specify the configuration file path in your workspace. Also update the output folder parameter "outputRoot" in the config file "example/input.1KG.GRCh37.txt" to specify a different output folder location.

To run CoCoRV using different inputs, you need to update "example/input.1KG.GRCh37.txt" and change case file specific parameters which are "caseBed", "caseVCFPrefix", "caseVCFSuffix", "caseSample" and output folder parameter which is "outputRoot".

The test data used here is also available to download from Amazon s3: s3://cocorv-1kg-grch37-data/
The processed gnomAD v2 exome data used here is also available to download from Amazon s3: s3://cocorv-resource-files/gnomADv2exome/

#### GRCh38 using gnmoAD v4 exome/genome data ####
To run the CoCoRV pipeline for GRCh38 using gnomAD v4 exome data or genome data, example run scripts "testGRCh38.gnomAD.v4exome_1KG_singularity.sh" and "testGRCh38.gnomAD.v4genome_1KG_singularity.sh" are given. These test scripts use the input configuration file given in here: "example/input.1KG.GRCh38.gnomAD.v4exomes.txt" and "input.1KG.GRCh38.gnomAD.v4genomes.txt". Here we used 23 test samples from 1000 Genomes Project (build GRCh38) as case data.

To run "testGRCh38.gnomAD.v4exome_1KG_singularity.sh" or "testGRCh38.gnomAD.v4genome_1KG_singularity.sh" script, update the "inputConfig" parameter to specify the configuration file path in your workspace. Also update the output folder parameter "outputRoot" in the config file in "example" to specify a different output folder location.

To run CoCoRV using different inputs, you need to update "example/input.1KG.GRCh38.gnomAD.v4exomes.txt" and change case file specific parameters which are "caseBed", "caseVCFPrefix", "caseVCFSuffix", "caseSample" and output folder parameter which is "outputRoot".

The test data used here is also available to download from Amazon s3: s3://cocorv-1kg-grch38-data/
The processed gnomAD v4 exome data used here is also available to download from Amazon s3: s3://cocorv-resource-files/gnomADv4exome/
The processed gnomAD v4 genome data used here is also available to download from Amazon s3: s3://cocorv-resource-files/gnomADv4genome/

### Contributors ###
* Saima Sultana Tithi
* Wenan Chen

### Contact ###
* Please contact Saima Sultana Tithi (saimasultana.tithi@stjude.org) or Wenan Chen (chen.wenan@mayo.edu) for any questions

### License ###
MIT license
