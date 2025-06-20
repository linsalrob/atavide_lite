# Metagenomics Analysis With Atavide

To address the limitations of monolithic metagenomic workflows, we developed `atavide_lite`, 
a modular, script-based pipeline that decouples computational processes to enhance robustness,
transparency, and adaptability across diverse high-performance computing (HPC) environments. 

This design allows users to inspect, rerun, or replace individual steps without affecting the 
overall pipeline logic or requiring global recomputation, features that are particularly 
advantageous when dealing with heterogeneous sequencing technologies, large sample sets, and unstable
job scheduling systems. 

The pipeline performs multiple steps that are standard for metagenomics processing. 

_Quality control and host read removal_. Read-based analyses begin with adapter trimming and quality
filtering using fastp, with parameters optimized for Illumina or Oxford Nanopore reads. 
Sequences are filtered to be longer than the minimum read length (100 bp) and we only
allow a single uncalled base. 

_Functional and taxonomic profiling_. The cleaned reads are converted to FASTA format using 
a compiled C utility. Taxonomic classification is performed using `mmseqs2 easy-taxonomy`, with 
UniRef50 as the default database. A summarisation pipeline, implemented in `snakemake`, adds taxonomic 
labels via `pytaxonkit`, allowing a more complete resolution of taxonomic lineages than native 
`mmseqs2` outputs. 

Taxonomy tables are generated at each rank for downstream comparative analyses.

Functional annotation leverages an SQLite database derived from the BV-BRC Subsystems framework.
Subsystem assignments are joined with the taxonomic tables to yield both raw and normalised counts, 
producing an integrated taxonomic-functional matrix.

_Assembly and metagenome-assembled genome (MAG) recovery_. _De novo_ assembly is carried
out for each sample independently using megahit, avoiding the memory constraints associated 
with cross-assemblies. Contigs from all samples are concatenated and binned using `vamb`, which 
employs a variational autoencoder for clustering based on tetranucleotide frequencies and 
abundance profiles. Read mapping for abundance estimation is performed using `minimap2`.

To validate MAG quality, we apply `checkm`, assessing completeness and contamination for bins derived
from both unsplit and split clusters. For stratified analyses 
(e.g. ecological or host-based subsets), we support grouped VAMB binning using a user-defined
TSV mapping of samples to groups, followed by independent execution of the concatenation, 
mapping, and binning steps per group.

