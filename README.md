# SuPreMo-Enformer

SuPreMo: Get sequences for predictive models

SuPreMo incorporates variants one at a time into the human reference genome and generates reference and alternate sequences for each perturbation under each provided augmentation parameter.
The sequences are accompanied by the relative position of the perturbation for each sequence.
The perturbation is by default centered in the sequence, if possible, unless the user provides a shifting parameter.
SuPreMo-Enformer: Get the disruption scores of gene expression and transcription factor binding using the Enformer model

SuPreMo-Enformer generates disruption scores by comparing contact frequency maps predicted from the reference and alternate sequences.
The maps are accompanied by the relative position of the perturbation and the start coordinate that the predicted maps correspond to.
The perturbation is by default centered in the map, if possible, unless the user provides a shifting parameter.
Optionally, the predicted maps or disruption score tracks along those maps can also be outputted.
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
