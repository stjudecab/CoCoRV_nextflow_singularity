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
Please look at Singularity documentation on how it install it on your system: https://docs.sylabs.io/guides/3.0/user-guide/installation.html

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
* *nextflow.config*: contains default configuration for running CoCoRV pipeline. Please update the "singularity.cacheDir" parameter inside "cluster_singularity" profile to specify a folder which can be used by Singularity to cache our pipeline's docker images (it can be any folder in your local cluster or local computer where you have full access).
* *example folder*: contains run specific configuration files. You need to create a configuration file like "input.1KG.GRCh38.gnomAD.v4.1exomes.txt" or "input.1KG.GRCh38.gnomAD.v4.1genomes.txt" given in here and update it with case VCF file path and output folder path for your run.

For detailed instructions on CoCoRV tool and it's parameters, please refer to the CoCoRV tool repository: https://bitbucket.org/Wenan/cocorv/src/master/

#### GRCh37 using gnmoAD v2 exome data ####
If you have GRCh37 data, you need to use this version. To run the CoCoRV pipeline for GRCh37 using gnomAD v2 exome data, an example run script "testGRCh37.gnomAD.v2exome_1KG_singularity.sh" is given. This test script uses the input configuration file given in here: "example/input.1KG.GRCh37.txt". Here we used 25 test samples from 1000 Genomes Project (build GRCh37) as case data.

To run "testGRCh37.gnomAD.v2exome_1KG_singularity.sh" script, update the "inputConfig" parameter to specify the configuration file path in your workspace. Also update the output folder parameter "outputRoot" in the config file "example/input.1KG.GRCh37.txt" to specify a different output folder location.

To run CoCoRV using different inputs, you need to update "example/input.1KG.GRCh37.txt" and change case file specific parameters which are "caseBed", "caseVCFPrefix", "caseVCFSuffix", "caseSample" and output folder parameter which is "outputRoot".

The test data used here is also available to download from Amazon s3 "[https://cocorv-1kg-grch37-data.s3.amazonaws.com/](https://cocorv-1kg-grch37-data.s3.amazonaws.com/)".
The processed gnomAD v2 exome data used here is also available to download from Amazon s3 "[https://cocorv-resource-files.s3.amazonaws.com/gnomADv2exome/](https://cocorv-resource-files.s3.amazonaws.com/gnomADv2exome/)".

#### GRCh38 using gnmoAD v4 exome/genome data ####
To run the CoCoRV pipeline for GRCh38 using gnomAD v4.1 exome data or genome data, example run scripts "testGRCh38.gnomAD.v4.1exome_1KG_singularity.sh" and "testGRCh38.gnomAD.v4.1genome_1KG_singularity.sh" are given. These test scripts use the input configuration file given in here: "example/input.1KG.GRCh38.gnomAD.v4.1exomes.txt" and "input.1KG.GRCh38.gnomAD.v4.1genomes.txt". Here we used 23 test samples from 1000 Genomes Project (build GRCh38) as case data.

To run "testGRCh38.gnomAD.v4.1exome_1KG_singularity.sh" or "testGRCh38.gnomAD.v4.1genome_1KG_singularity.sh" script, update the "inputConfig" parameter to specify the configuration file path in your workspace. Also update the output folder parameter "outputRoot" in the config file in "example" to specify a different output folder location.

To run CoCoRV using different inputs, you need to update "example/input.1KG.GRCh38.gnomAD.v4.1exomes.txt" and change case file specific parameters which are "caseBed", "caseVCFPrefix", "caseVCFSuffix", "caseSample" and output folder parameter which is "outputRoot".

The test data used here is also available to download from Amazon s3 "[https://cocorv-1kg-grch38-data.s3.amazonaws.com/](https://cocorv-1kg-grch38-data.s3.amazonaws.com/)".
The processed gnomAD v4 exome data used here is also available to download from Amazon s3 "[https://cocorv-resource-files.s3.amazonaws.com/gnomADv4.1exome/](https://cocorv-resource-files.s3.amazonaws.com/gnomADv4.1exome/)".
The processed gnomAD v4 genome data used here is also available to download from Amazon s3 "[https://cocorv-resource-files.s3.amazonaws.com/gnomADv4.1genome/](https://cocorv-resource-files.s3.amazonaws.com/gnomADv4.1genome/)".

#### Download gnomAD data using AWS CLI ####
As this is a huge dataset, it is better to use Amazon AWS command line tool aws-cli to download the data. 

Here is how you can install aws-cli:
https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html

After installing this, you can use "aws s3" command to list any s3 bucket folder, or download any folder or files from s3.
https://docs.aws.amazon.com/cli/latest/reference/s3/

Here are the s3 bucket paths of the gnomAD data:
s3://cocorv-resource-files/gnomADv4.1exome/
s3://cocorv-resource-files/gnomADv4.1genome/

To download the data, you need to run commands like this:
```bash
cd /local-dir-path-where-you-want-download/
aws s3 cp s3://cocorv-resource-files/gnomADv4.1exome/ . --recursive
```
You can check all resource files for CoCoRV using this command:
```bash
aws s3 ls s3://cocorv-resource-files/
```

#### TODO ####
* Move local path specific profiles to the input configuration file

### Contributors ###
* Saima Sultana Tithi
* Wenan Chen

### Contact ###
* Please contact Saima Sultana Tithi (saimasultana.tithi@stjude.org) or Wenan Chen (chen.wenan@mayo.edu) for any questions

### License ###
MIT license
