# A nextflow based pipeline for CoCoRV based analysis #

### What is this repository for? ###
This repository includes the nextflow based pipeline for applying CoCoRV to sequencing data, e.g., whole exome sequencing (WES), or whole genome sequencing (WGS), using gnomAD data as control. It supports gnomAD v2 GRCh37 exome data, gnomAD v3 GRCh38 whole genome data, and gnomAD v4 exome data as control.    

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
```bash
nextflow run CoCoRVPipeline.nf
```
It will print the help message.

### Examples ###
The CoCoRV nextflow pipeline script, nextflow configuration file, and example run configuration files are provided in this repository.

* *CoCoRVPipeline.nf*: main nextflow script for CoCoRV pipeline.
* *nextflow.config*: contains default configuration for running CoCoRV pipeline.
* *example folder*: contains run specific configuration files. You need to create a configuration file like "input.GRCh38.gnomAD.v4exomes.txt" given in here and update it with case VCF file path and output folder path for your run.

#### GRCh38 using gnmoAD v3 genome data ####
To run the CoCoRV pipeline for GRCh38 using gnomAD v4 exome data, an example run script "test_GRCh38.gnomAD.v4exomes_singularity.sh" is given. This test script uses the input configuration file given in here: "example/input.GRCh38.gnomAD.v4exomes.txt". Here we used 407 ALS sporadic case samples (build GRCh38) as case data.

To run "test_GRCh38.gnomAD.v4exomes_singularity.sh" script, update the "inputConfig" parameter to specify the configuration file path in your workspace. Also update the output folder parameter "outputRoot" in the config file "example/input.GRCh38.gnomAD.v4exomes.txt" to specify a different output folder location.

To run CoCoRV using different inputs, you need to update "example/input.GRCh38.gnomAD.v4exomes.txt" and change case file specific parameters which are "caseBed", "caseVCFPrefix", "caseVCFSuffix", "caseSample" and output folder parameter which is "outputRoot".

### Contributors ###
* Saima Sultana Tithi
* Wenan Chen

### Contact ###
* Please contact Saima Sultana Tithi (saimasultana.tithi@stjude.org) or Wenan Chen (chen.wenan@mayo.edu) for any questions

### License ###
MIT license
