# SuPreMo-Enformer

SuPreMo: Get sequences for predictive models

SuPreMo incorporates variants one at a time into the human reference genome and generates reference and alternate sequences for each perturbation under each provided augmentation parameter.
The sequences are accompanied by the relative position of the perturbation for each sequence.
The perturbation is by default centered in the sequence, if possible, unless the user provides a shifting parameter.
SuPreMo-Enformer: Get the disruption scores of gene expression and transcription factor binding using the Enformer model

SuPreMo-Enformer generates disruption scores by comparing RNA-seq, ChIP-seq, CAGE-seq and DNASE tracks predicted from the reference and alternate sequences.

# Supremo Enformer

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](link-to-build-status)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](CHANGELOG.md)

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#Usage)
- [Data Requirements and Preprocessing](#data-requirements-and-preprocessing)
- [Output Format](#output_format)
- [Examples and Tutorials](#examples-and-tutorials)
- [Citation and References](#citation-and-references)
- [License](#license)
- [Contact and Support](#contact-and-support)
- [Acknowledgments](#acknowledgments)

## Introduction
**Supremo Enformer** is a deep learning framework designed to predict genomic regulatory functions directly from DNA sequence data. Leveraging a transformer-based architecture, Supremo Enformer captures long-range interactions and integrates diverse regulatory signals, such as enhancer activity, transcription factor binding, and gene expression. This tool is ideal for large-scale genomic analyses and detailed investigations into the regulatory mechanisms that drive gene expression.

## Features
- **Transformer-based architecture:** Captures long-range genomic dependencies.
- **Integrated predictions:** Simultaneously predicts multiple regulatory elements.
- **Scalability:** Optimized for large datasets and high-throughput analysis.
- **Enhanced interpretability:** Provides insights into genomic regulatory mechanisms.

## Installation
### Create conda environment
conda env create -f supremo_enformer.yml \
conda activate supremo_enformer_env

### Installation Steps
1. **Clone the repository:**
   ```bash
   git clone https://github.com/May-BG/SuPreMo-Enformer.git
   cd supremo_enformer


## Usage
```shell
python streamlined_SuPreMo.py \
    "$data_dir$vcf_file" \
    --dir "$data_dir/results" \
    --file "skin_melanoma_DEL_part_$1.out" \
    --get_Enformer_scores \
--selected_tracks $(awk -F'\t' '{print $NF}' target_new.txt | sort | uniq)
```

### Command-Line Arguments
--file: name of outputs. \
--dir: Directory to save results. \
--selected_tracks: (Optional) Path to the selected track names from Enformer. \
--get_Enformer_scores: required.

## data-requirements-and-preprocessing
- VCF file
    * following [vcf 4.1/4.2 specifications](https://samtools.github.io/hts-specs/VCFv4.1.pdf)
- TXT file
    * Columns required for simple variants: CHROM, POS, REF, ALT
    * Columns required for structural variants: CHROM, POS, REF, ALT, END, SVTYPE (SV type), SVLEN (SV length)
 
```shell
python /pollard/data/projects/xzhang/tcga/SuPreMo/enformer/SuPreMo-enformer/manual_load_enformer_scripts/streamlined_SuPreMo.py VARIANT_FILE --dir OUTPUT_DIRECTORY --file OUTPUT_NAME --get_Enformer_scores --cell_type CELL_TYPE_FROM_TARGET_FILE --protocol ATAC/CHIP/CAGE/DNASE
```
Note: Replace the words in ALL CAPS with custom values:
- VARIANT_FILE: Path to your variant file (e.g., /pollard/data/projects/xzhang/tcga/del_100.txt).
- OUTPUT_DIRECTORY: Directory for saving output (e.g., /pollard/data/projects/xzhang/tcga).
- OUTPUT_NAME: Name for the output file (e.g., del_100_new.out).
--get_Enformer_scores: Use this flag as-is if you want Enformer scores.
--cell_type CELL_TYPE_FROM_TARGET_FILE: Replace with the cell type from targets_human.txt (e.g., "K562").
--protocol ATAC/CHIP/CAGE/DNASE: Specify the protocol, such as ATAC, CHIP, CAGE, or DNASE.
* Conditional Argument:
If --protocol is CHIP, add the --binding_factor argument (e.g., --binding_factor CTCF).
If --protocol is not CHIP, do not include --binding_factor.


### example:
```shell
python /pollard/data/projects/xzhang/tcga/SuPreMo/enformer/SuPreMo-enformer/manual_load_enformer_scripts/streamlined_SuPreMo.py /pollard/data/projects/xzhang/tcga/del_1.txt --dir /pollard/data/projects/xzhang/tcga --file del_1_dnase_colon.out --get_Enformer_scores --cell_type "colon epithelial cell line" --protocol DNASE

```

In this example, --binding_factor is omitted because --protocol is DNASE.
