# A nextflow based pipeline for CoCoRV based analysis #

### What is this repository for? ###
This repository includes the nextflow based pipeline for applying CoCoRV to sequencing data, e.g., whole exome sequencing (WES), or whole genome sequencing (WGS), using gnomAD data as control. It supports both gnomAD v2 GRCh37 exome data and gnomAD v3 GRCh38 whole genome data as control.    

### Installation ###
#### Micromamba ####
Install Micromamba using the following instructions:
```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
```
It will download the micromamba binary in "bin" folder.

Update your home .bashrc using the following command, so that you can run the micromamba command from anywhere.
```bash
./bin/micromamba shell init -s bash -p ~/micromamba  # this writes to your home .bashrc file and creates a micromamba folder in your home directory which will contain all of your conda environments created by micromamba
```
Restart your terminal and check your micromamba installtion by running micromamba command.
```bash
micromamba
```

If you need to install micromamba in macOS or Windows, please check Micromamba installation page: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

#### Create Micromamba environment for gnomAD v3 GRCh38 ####
Download the git repository "stjudecab/CoCoRV_nextflow_conda". Go to "conda-envs" folder and create conda environment "cocorv-env-GRCh38":
```bash
git clone https://github.com/stjudecab/CoCoRV_nextflow_conda.git
cd CoCoRV_nextflow_conda/conda-envs
micromamba create -f cocorv-env-GRCh38.yaml
```
It will create a conda environment with name "cocorv-env-GRCh38" with all necessary python packages, gnomAD packages, and R packages. 
Activate the "cocorv-env-GRCh38" environment and install custom CoCoRV R package in there.
```bash
micromamba activate cocorv-env-GRCh38
cd /home/stithi/CoCoRV_pipeline_bitbucket/cocorv
Rscript build.R
micromamba deactivate
```

#### Create Micromamba environment for gnomAD v2 GRCh37 ####
Download the git repository "stjudecab/CoCoRV_nextflow_conda". Go to "conda-envs" folder and create conda environment "cocorv-env-GRCh37":
```bash
git clone https://github.com/stjudecab/CoCoRV_nextflow_conda.git
cd CoCoRV_nextflow_conda/conda-envs
micromamba create -f cocorv-env-GRCh37.yaml
```
It will create a conda environment with name "cocorv-env-GRCh37" with all necessary python packages, gnomAD packages, and R packages. 
Activate the "cocorv-env-GRCh37" environment and install custom CoCoRV R package in there.
```bash
micromamba activate cocorv-env-GRCh37
cd /home/stithi/CoCoRV_pipeline_bitbucket/cocorv
Rscript build.R
micromamba deactivate
```

#### Nextflow ####
Load nextflow using the following commands:
```bash
module load nextflow/23.04.1
module load java/openjdk-11
```
If you want to install nextflow, please check https://www.nextflow.io/.

#### Run CoCoRV pipeline ####
```bash
nextflow run CoCoRVPipeline.nf
```
It will print the help message.

### Examples ###
The CoCoRV nextflow pipeline script, nextflow configuration file, and example run configuration files are provided in this repository.

* *CoCoRVPipeline.nf*: main nextflow script for CoCoRV pipeline.
* *nextflow.config*: contains default configuration for running CoCoRV pipeline, please update the "process.conda" parameter to specify your home directory's micromamba environment folder containing "cocorv-env-GRCh38" or "cocorv-env-GRCh37" environment.
* *example folder*: contains run specific configuration files. You need to create a configuration file like "input.1KG.GRCh38.txt" given in here and update it with case VCF file path and output folder path for your run.
* *conda-envs folder*: contains YAML files for creating conda environment.

#### GRCh38 using gnmoAD v3 genome data ####
To run the CoCoRV pipeline for GRCh38 using gnomAD v3 whole genome data, an example run script "testGRCh38_1KG_conda.sh" is given. This test script uses the input configuration file given in here: "example/input.1KG.GRCh38.txt". Here we used 23 samples from 1000 Genomes Project whole genome data (build GRCh38) as case data.

To run "testGRCh38_1KG_conda.sh" script, update the "inputConfig" parameter to specify the configuration file path in your workspace. Also update the output folder parameter "outputRoot" in the config file "example/input.1KG.GRCh38.txt" to specify a different output folder location. Also update "process.conda" parameter in "nextflow.config" file's *conda_GRCh38* profile section to specify your home directory's micromamba environment folder containing "cocorv-env-GRCh38" environment (usually YOUR_HOME_DIRECTORY/micromamba/envs/cocorv-env-GRCh38).

To run CoCoRV using different inputs, you need to update "example/input.1KG.GRCh38.txt" and change case file specific parameters which are "caseBed", "caseVCFPrefix", "caseVCFSuffix", "caseSample" and output folder parameter which is "outputRoot".

#### GRCh37 using gnmoAD v2 exome data ####
To run the CoCoRV pipeline for GRCh37 using gnomAD v2 exome data, an example run script "testGRCh37_1KG_conda.sh" is given. This test script uses the input configuration file given in here: "example/input.1KG.GRCh37.txt". Here we used 25 samples from 1000 Genomes Project whole exome data (build GRCh37) as case data.

To run "testGRCh37_1KG_conda.sh" script, update the "inputConfig" parameter to specify the configuration file path in your workspace. Also update the output folder parameter "outputRoot" in the config file "example/input.1KG.GRCh37.txt" to specify a different output folder location. Also update "process.conda" parameter in "nextflow.config" file's *conda_GRCh37* profile section to specify your home directory's micromamba environment folder containing "cocorv-env-GRCh37" environment (usually YOUR_HOME_DIRECTORY/micromamba/envs/cocorv-env-GRCh37).

To run CoCoRV using different inputs, you need to update "example/input.1KG.GRCh37.txt" and change case file specific parameters which are "caseBed", "caseVCFPrefix", "caseVCFSuffix", "caseSample" and output folder parameter which is "outputRoot".

### Contributors ###
* Saima Sultana Tithi
* Wenan Chen

### Contact ###
* Please contact Saima Sultana Tithi (saimasultana.tithi@stjude.org) or Wenan Chen (chen.wenan@mayo.edu) for any questions

### License ###
MIT license
