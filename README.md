# Rovaca

## Introduction

Rovaca is a tool for detecting SNPs (Single Nucleotide Polymorphisms) and INDELs (Insertions/Deletions) from DNA sequencing data in a single sample. It is a probabilistic model-based variant detection algorithm designed to accurately identify variants in samples.

The working principles of Rovaca are as follows:

1. **Identifying active regions**: The program determines which regions of the genome (active regions) to process based on the presence of evidence for variants.
2. **Determining haplotypes through assembly of active regions**: Rovaca uses a local de-novo assembly method to split sequencing data into smaller fragments (haplotypes). These fragments are DNA segments caused by potential variants. The program then uses the Smith-Waterman algorithm to realign each haplotype with the reference haplotype to identify potential variant sites.
3. **Pairwise alignment using PairHMM**: For each active region, the program performs pairwise alignment for each read of each haplotype using the PairHMM algorithm. This produces a likelihood matrix of haplotypes given the read data. These likelihoods are then marginalized to obtain the likelihood of each allele at each potential variant site given the read data.
4. **Variant calling using Bayesian inference**: For each potential variant site, the program applies Bayes' rule, using the likelihood of alleles given the read data to calculate the likelihood of each genotype for each sample given the read data for that sample. The most likely genotype is then assigned to the sample.

Rovaca offers high sensitivity and specificity, capable of accurately detecting variants in complex genomic regions. It can also process multi-sample data, performing joint variant detection to improve accuracy.

## Prerequisites

- **System requirements**: Intel or AMD x86 system with at least **AVX2** support.
- **For best performance**: AVX512 is recommended.
- **Checking compatibility**: Run the following command to check AVX support:

  ```sh
  lscpu | grep avx
  ```

## Installation

We strongly recommend users to **download the released version** for ease of use. However, if you prefer to compile from the source code, follow these steps:

1. Ensure you have **libboost** version **>=1.69.0** installed.
2. Clone the repository and navigate to the source directory.
3. Use the provided script to build the software:

   ```sh
   ./build_dev.sh
   ```

## Usage

```sh
# Basic usage with required parameters
rovaca <tool> <options>
```

Currently, only HaplotypeCaller is supported as a tool; other tools may be integrated as needed in the future.

### Parameters

```sh
-H, --help                                    Display help information
-V, --version                                 Display version information
-I, --input <file1> <file2> ...               Input file path(s), multiple files can be specified
-O, --output <file>                           Output file path, required parameter
-R, --reference <file>                        Reference file path, required parameter
-L, --interval <file>                         Interval file path
-P, --interval-padding <value>                Size to expand intervals on both sides, non-negative integer, default: 0
-D, --max-reads-depth <value>                 Maximum number of reads at each alignment start position, non-negative integer, default: 50
-Q, --base-quality-score-threshold <value>    Minimum base quality score threshold, [6, 127], default: 18
-G, --gvcf-gq-bands <value1> <value2> ...     Specify GQ boundaries for merging non-variant sites in GVCF mode, [1, 100], default: [1, 2, 3,... 60, 70, 80, 90, 99]
--nthreads <value>                            Number of threads to start, [1, 128], default: 30
--pcr-indel-model <value>                     PCR indel model, default: CONSERVATIVE, options: {NONE, HOSTILE, CONSERVATIVE, AGGRESSIVE}
--emit-ref-confidence <value>                 Reference confidence model, default: NONE, options: {NONE, GVCF}
--nstreampool <value>                         Iostream pool size, [1, 20], default: 10, usually does not need to be specified
--inspect-reads                               Strictly verify each read in the input BAM file, default: false
--bqsr-recal-table <file>                     Specify the recalibration.table file for Base Quality Score Recalibration (BQSR)
--compression-level <value>                   Compression level of the output file, only effective for gz files, [0-9], default: 6
```

## Examples

Here are some example commands for running Rovaca:

```sh
# Basic usage with required parameters
rovaca HaplotypeCaller -I sample.bam -O output.vcf -R reference.fasta

# Running with interval padding and multiple threads
rovaca HaplotypeCaller -I sample.bam -O output.vcf -R reference.fasta -P 100 --nthreads 16

# Generating GVCF output
rovaca HaplotypeCaller -I sample.bam -O output.g.vcf -R reference.fasta --emit-ref-confidence GVCF
```

## Additional Information

For further details and updates, please refer to the official documentation or reach out to the development team.
