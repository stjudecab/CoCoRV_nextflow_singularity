# nf-core/rarevariantburden: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rarevariantburden/usage](https://nf-co.re/rarevariantburden/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-core/rarevariantburden** is a bioinformatics pipeline that performs consistent summary count based rare variant burden test, which is useful when we only have sequenced cases/patients data, no matched control data, here we provided pre-processed and annotated public summary count data, such as gnomAD data, which can be used for rare variant burden test and can be used to identify disease-predisposition genes present in the case study.

## Prerequisites

1. Install Nextflow (>=23.10.0) using the instructions [here.](https://nextflow.io/docs/latest/getstarted.html#installation)
2. Install one of the following technologies for full pipeline reproducibility: Docker, Singularity, Podman, Shifter or Charliecloud.

## Run nf-core/rarevariantburden with test data

Before running the pipeline with your data, we recommend running it with the test dataset available [here](https://github.com/nf-core/test-datasets/tree/rarevariantburden). You do not need to download the data as the pipeline is configured to fetch that data automatically for you when you use the test profile.

Run the following command, where YOURPROFILE is the package manager you installed on your machine. For example, `-profile test,docker` or `-profile test,singularity`:

```bash
nextflow run nf-core/rarevariantburden \
    -profile test,<YOURPROFILE> \
    --outdir <OUTDIR>
```

> Check [nf-core/configs](https://github.com/nf-core/configs/tree/master/conf) to see if a custom config file to run nf-core pipelines already exists for your institute. If so, you can simply use `-profile test,<institute>` in your command. This enables the appropriate package manager and sets the appropriate execution settings for your machine.
> NB: The order of profiles is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

Running the command creates the following files in your working directory:

```
work                # Directory containing the Nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other Nextflow hidden files, like history of pipeline logs.
```

## Running the pipeline with your data
Running the pipeline involves three steps:

1. Prepare joint called VCF file:
  * First, prepare the joint called and VQSR applied VCF file from your case study. You can use [nf-core/sarek](https://nf-co.re/sarek/) pipeline's [GATK joint calling module](https://nf-co.re/sarek/3.5.1/docs/output/#gatk-joint-germline-variant-calling) to prepare a joint called and VQSR applied VCF file from your sample VCF files. You also need to prepare a text file containing sample IDs, one sample ID per line.

2. Download control data:
  * For control data, you need to download the control data from our Amazon AWS s3 bucket. We provide 3 different control datasets, For build GRCH37, we have gnomAD v2 exome data, for build GRCh38, we have gnomAD v4.1 exome and gnomAD v4.1 genome data as controls.
  * As the control data is a huge dataset, it is better to use Amazon AWS command line tool [aws-cli](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) to download the data. After installing this, you can use "aws s3" command to list any s3 bucket folder, or download any folder or files from s3. You will find the s3 commands list [in here](https://docs.aws.amazon.com/cli/latest/reference/s3/).
  * Here are the s3 bucket paths of the 3 gnomAD control datasets:
    * s3://cocorv-resource-files/gnomADv2exome/
    * s3://cocorv-resource-files/gnomADv4.1exome/
    * s3://cocorv-resource-files/gnomADv4.1genome/
  * To download the data, you need to run following command:
  ```bash
  cd /local-dir-path-where-you-want-download/
  aws s3 cp s3://cocorv-resource-files/gnomADv2exome/ . --recursive
  ```
  * You can check all resource files for our pipeline using this command:
  ```bash
  aws s3 ls s3://cocorv-resource-files/
  ```

3. Download the annotation tool resources:
  * You also need to download the annovar and VEP resource folders for running Annovar and VEP annotation.
  * Here are the s3 bucket paths of the annotation tool datasets:
    * s3://cocorv-resource-files/annovarFolder/
    * s3://cocorv-resource-files/vepFolder/

Now, you can run the pipeline using the following command:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/rarevariantburden \
   -profile <docker/singularity/.../institute> \
   --caseJointVCF <jointVCF.vcf.gz> \
   --caseSample <sampleList.txt> \
   --controlDataFolder <controldataFolder> \
   --annovarFoler <annovarFolder> \
   --vepFolder <vepFolder> \
   --build <GRCh37/GRCh38> \
   --gnomADVersion <v2exome/v4exome/v4genome> \
   --outdir <OUTDIR>
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/rarevariantburden -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
caseJointVCF: 'jointVCF.vcf.gz'
caseSample: 'sampleList.txt'
controlDataFolder: 'controldataFolder'
annovarFolder: 'annovarFolder'
vepFolder: 'vepFolder'
build: 'GRCh37'
gnomADVersion: 'v2exome'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/rarevariantburden
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rarevariantburden releases page](https://github.com/nf-core/rarevariantburden/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
