# QSPipeline

QSPipeline is a pipeline to detect viral quasispecies in NGS data. It is based in R and can be run via RStudio or via command-line script. It is compatible with Linux, Mac OS and Windows.

## Table of contents

1. [The QSPipeline workflow](#workflow)
2. [Dependencies](#dependencies)
3. [Installation](#installation)
4. [Usage](#usage)

## The QSPipeline workflow <a name="workflow"></a>

FIGURE

The image above shows the workflow the user should follow to obtain the viral quasispecies present in its samples:

  1. Check metadata: to verify primer descriptors and files for pools in run.
<!-- -->
  2. Quality assesment: to analyze samples' quality.
<!-- -->
  3. QSPipeline: to obtain quasispecies present in the samples after applying quality, similarity and abundance filters.

## Dependencies <a name="dependencies"></a>

To run QSPipeline, several R libraries are required:

  * [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
  * [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html)
  * [data.table](https://cran.r-project.org/web/packages/data.table/)
  * [stringr](https://cran.r-project.org/web/packages/stringr/)
  * [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/)
  * [optparse](https://cran.r-project.org/web/packages/optparse/)

FLASH of Illumina is required too, but it is provided in this repository.

## Installation <a name="installation"></a>

Clone the QSPipeline repository in directory of your choice:

```
git clone https://github.com/bsorianos/QSPipeline.git
```

Or download it directly from GitHub QSPipeline page.

This repository includes:

  * a folder called FLASH-1.2.11, which includes FLASH of Illumina.<br />
  * a folder called R, which includes several scripts needed to run the pipeline.
  * the data folder containing three files: samples.csv, primers.csv and AmpliconRefSeqs.fna.
  * the export folder containing some pipeline results.
  * a folder called filt, which includes filtered haplotypes.
  * the flash folder containing the files resulting from FLASH.
  * the flashFilt folder including filtered files by quality.
  * a folder called join, which contains the results coming from intersecting haplotypes between FW and RV files.
  * the nt.01 folder containing final haplotypes per region with an abundance filter of 0.1% by default.
  * the nt folder containing final haplotypes per region with an abundance filter of 0.5% by default.
  * a folder called reports including every report generated in each step of the pipeline.
  * the results folder containing a summary file of each sample and a coverage boxplot by amplicon.
  * the run folder containing the FASTQ files that are going to be analyzed.
  * the trim folder containing fasta files after adapter trimming.
  * the 00_CheckMetadata-RAVs-v1.17.R file.
  * the 01_MiSeq_RAV_QA_Pipeline-v2.2.R file.
  * the 02_QSPipeline-v1.06.R file.
  * this README.md file.
  * the LICENSE file.
  * the default_pars.R file, which contains the default parameters used by this pipeline.

## Usage <a name="usage"></a>

Follow these steps to execute the complete pipeline:

1. Run the 00_CheckMetadata-RAVs-v1.17.R file to verify all primer descriptors and files for pools.
<!-- -->
2. Run the 01_MiSeq_RAV_QA_Pipline-v2.2.R file to asses the quality of the R1 and R2 fastq files and filter reads by quality.
<!-- -->
3. Run the 02_QSPipeline-v1.06.R file to detect the quasispecies present in the samples. To configure this pipeline, it has the following options:

    * -p/--pmm_mx: Maximum number of mismatches in the specific primer.
    * -l/--min_len: Minimum length to consider a sequence.
    * -m/--min_reads: Minimum number of reads by sequences after repair.
    * -n/--max_ns: Maximum number of admissible Ns.
    * -d/--max_diffs: Maximum number of tolerated differences between the haplotype and its reference.
    * -g/--max_gaps: Maximum number of admissible gaps.
    * -r/--ref_type: Reference type to filter reads. Choose between generic or consensus.
    * -M/--method: Sum or Intersect. Sum takes the sum of the common haplotypes reads as a distribution, while Intersect takes the intersection as a distribution.
    * -t/--var_thr: Accept variant with an abundance above this value.
    * -a/--ab_thr: Second abundance filter for haplotypes.
    * -s/--min_size: Value between 0 and 1, which multiplies sequence length, to select minimum length of sequences.

To run these R scripts, you can execute them via terminal (for example, Rscript 02_QSPipeline-v1.06.R --max_gaps 50) or you can run them via RStudio (for example, system('Rscript 02_QSPipeline-v1.06.R --max_gaps 50').
