<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-rarevariantburden_logo_dark.png">
    <img alt="nf-core/rarevariantburden" src="docs/images/nf-core-rarevariantburden_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/rarevariantburden/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/rarevariantburden/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/rarevariantburden/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/rarevariantburden/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rarevariantburden/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/rarevariantburden)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rarevariantburden-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rarevariantburden)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

#### TOC

- [Introduction](#introduction)
- [Pipeline summary](#pipeline-summary)
- [Usage](#usage)
- [Pipeline output](#pipeline-output)
- [Credits](#credits)
- [Contributions and Support](#contributions-and-support)
- [Citations](#citations)

## Introduction

**nf-core/rarevariantburden** is a bioinformatics pipeline that performs consistent summary count based rare variant burden test, which is useful when we only have sequenced cases/patients data, no matched control data, here we provided pre-processed and annotated public summary count data, such as gnomAD data, which can be used for rare variant burden test and can be used to identify disease-predisposition genes present in the case study.

Some key features of our pipeline:

- Consistent filtering is applied to make sure the same set of high quality variants are used.
- It can stratify cases into different ethnicity groups, and perform stratified analysis with group-matched control summary counts.
- For recessive models, it can exclude double heterozygous due to high linkage disequilibrium in populations.
- Also provides accurate inflation factor estimate, QQ plot, and powerful FDR control for discrete count data, whose p-value distribution under the null is far from the uniform distribution when the alleles are very rare.
- It supports gnomAD v2 exome (GRCh37) data, and gnomAD v4.1 exome (GRCh38) data, and gnomAD v4.1 genome (GRCh38) data as control.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Pipeline summary

<picture align="center">
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/cocorv_subway.png">
    <img alt="nf-core/rarevariantburden workflow diagram" src="docs/images/cocorv_subway.png">
</picture>

1. Split the case joint called and VQSR applied VCF files chromosomewise (Using [BCFtools](https://samtools.github.io/bcftools/bcftools.html))
2. Normalize and QC the splitted case VCF files (Using [BCFtools](https://samtools.github.io/bcftools/bcftools.html))
3. Annotate normalized and QC'd VCF files with [Annovar](https://annovar.openbioinformatics.org/en/latest/) and [VEP](https://www.ensembl.org/vep) (VEP annotation is optional)
4. Convert the normalized and annotated VCF files to GDS format, which is easier to process in R (Using R seqarray)
5. Predict the ethnicity of the case samples (Using gnomAD random forest classifier)
6. Perform assiciation test for each VCF file using our [CoCoRV](https://bitbucket.org/Wenan/cocorv/) (Consistent summary Count based Rare Variant burden test) R package
7. Merge association test results
8. Calculate false positive rate (FDR) from merged results, plot QQ plot and lambda value using different R libraries
9. For top K genes, generate the list of samples and associated variants along with the annotations for the variants, this list will help the users to further check the top genes and their variants

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

First, prepare the joint called and VQSR applied VCF file from your case study. You can use [nf-core/sarek](https://nf-co.re/sarek/) pipeline's [GATK joint calling module](https://nf-co.re/sarek/3.5.1/docs/output/#gatk-joint-germline-variant-calling) to prepare a joint called and VQSR applied VCF file from your sample VCF files. You also need to prepare a text file containing sample IDs, one sample ID per line.

For control data, you need to download the control data from our Amazon AWS s3 bucket. We provide 3 different control datasets, For build GRCH37, we have gnomADv2exome data, for build GRCh38, we have gnomADv4.1exome and gnomADv4.1genome data as controls.

As the control data is a huge dataset, it is better to use Amazon AWS command line tool [aws-cli](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) to download the data.

After installing this, you can use "aws s3" command to list any s3 bucket folder, or download any folder or files from s3. You will find the s3 commands list [in here](https://docs.aws.amazon.com/cli/latest/reference/s3/).

Here are the s3 bucket paths of the 3 gnomAD control datasets:

- s3://cocorv-resource-files/gnomADv2exome/
- s3://cocorv-resource-files/gnomADv4.1exome/
- s3://cocorv-resource-files/gnomADv4.1genome/

To download the data, you need to run following command:

```bash
cd /local-dir-path-where-you-want-download/
aws s3 cp s3://cocorv-resource-files/gnomADv2exome/ . --recursive
```

You can check all resource files for our pipeline using this command:

```bash
aws s3 ls s3://cocorv-resource-files/
```

You also need to download the annovar and VEP resource folders for running Annovar and VEP annotation.

Here are the s3 bucket paths of the annotation tool datasets:

- s3://cocorv-resource-files/annovarFolder/
- s3://cocorv-resource-files/vepFolder/

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

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/rarevariantburden/usage) and the [parameter documentation](https://nf-co.re/rarevariantburden/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rarevariantburden/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/rarevariantburden/output).

## Credits

nf-core/rarevariantburden is written by Saima Sultana Tithi (saimasultana.tithi@stjude.org) and Wenan Chen (chen.wenan@mayo.edu).

<!-- TODO nf-core:
We thank the following people for their extensive assistance in the development of this pipeline:
If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#rarevariantburden` channel](https://nfcore.slack.com/channels/rarevariantburden) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/rarevariantburden for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

To learn more about the original CoCoRV tool, please look at our paper published in _Nature Communications_ [Pubmed link](https://pubmed.ncbi.nlm.nih.gov/35545612/).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
