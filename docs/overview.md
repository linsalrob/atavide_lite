# Overview

## What problem does atavide lite solve?

Metagenomics datasets routinely contain hundreds or thousands of samples. A typical analysis performs quality filtering, host screening, taxonomic and functional annotation, assembly, read mapping, binning, and MAG quality assessment. A failure in any sample or any stage can interrupt a monolithic workflow and make it difficult to see which command, input, dependency, or resource limit caused the problem.

HPC portability is a second problem. Two clusters using the same scheduler can still differ in partition names, maximum wall times, array limits, memory accounting, local and shared scratch, filesystem purge policies, software modules, GPU directives, and project accounting. A resource request optimised for one system may be invalid or inefficient on another.

`atavide_lite` responds by exposing the analysis as a sequence of small batch scripts. Each stage:

- has a clear input and output;
- can be inspected before submission;
- writes scheduler logs that identify the failing sample or job;
- can be rerun without repeating successful upstream work;
- can use scheduler dependencies to preserve ordering; and
- can request resources appropriate to both the task and cluster.

## Supported analyses

The main profiles cover two input models:

- **Paired-end short reads**, normally named with `_R1` and `_R2`, from platforms such as Illumina or MGI.
- **Single-end or long reads**, represented by one FASTQ file per sample, including Oxford Nanopore data.

The workflow combines two complementary branches:

- **Read-based analysis** retains per-read abundance information and produces taxonomic and BV-BRC Subsystems summaries.
- **Assembly-based analysis** creates contigs with MEGAHIT, maps reads back to contigs, groups contigs into MAGs with VAMB, and can assess bins with CheckM.

Host removal is optional. The “host” can be any reference sequence or multi-FASTA collection that should be separated from downstream analysis; it is not restricted to human data.

## What the project provides

The repository contains:

- cluster-specific Slurm profiles for Pawsey Setonix and Flinders Deepthought;
- PBS scripts for NCI Gadi;
- paired-read and single/long-read variants;
- adapter sequences;
- Conda/Mamba environment definitions;
- small compiled and Python utilities;
- Snakemake subworkflows for selected summarisation tasks; and
- documentation and profile-specific run books.

The profile scripts are intended as readable operational recipes. They are not a promise that every default resource request is optimal for every dataset. Users should start small, inspect scheduler usage, and tune time, memory, CPU, and storage requests within their local policy.

## Scope and responsibilities

`atavide_lite` coordinates established bioinformatics tools; it does not replace understanding their outputs or limitations. Users remain responsible for:

- ethical and authorised handling of sequence and host data;
- selecting appropriate references and databases;
- checking quality-control and host-removal results;
- validating taxonomic, functional, assembly, and MAG interpretations;
- complying with HPC allocation and storage policies; and
- citing the project and the tools used in their analysis.
