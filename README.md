[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.sdsu.edu/research)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/atavide_lite)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15356766.svg)](https://doi.org/10.5281/zenodo.15356766)

# atavide lite

A simpler version of atavide that relies only on slurm or PBS scripts. Some of the settings maybe specific for our compute resources

[atavide](https://github.com/linsalrob/atavide) is a full-fledged metagenome processing pipeline. This is a simpler version that uses either slurm scripts or PBS scripts to accomplish the same things. In this initial release, we avoided snakemake to work through the steps one at a time, although that will likely come soon.

# Pipeline steps

1. Run fastp to trim Illumina adapters
2. Use `minimap2` and `samtools` to filter out host and not host reads. Currently host reads are ignored. Host can be human or anything else.
3. Use `mmseqs easy-taxonomy` to compare not host reads to UniRef. By default we use UniRef50 but you could use any version.
4. Create a taxonomic summary for each sample and make a single `.tsv` file.
5. Connect in the subsystems from BV-BRC, and make a table that includes subsystems and taxonomic information
6. Create a subsystems taxonomy for the data


# different versions

## Paired vs Single End

In our modern processing, we have:
 - Paired end reads from MGI or Illumina sequencing. Those files usually end `_R1.fastq.gz` and `_R2.fastq.gz`, and the code looks for those. 
 - Single end reads from ONT sequencing. These files end `.fastq.gz`

We have two versions of the pipeline that work with _either_ paired end or single end. 

**However:** If you download some sequences from SRA, ENA, or DDBJ they may have paired end reads that end `_1.fastq.gz` and `_2.fastq.gz` in which case you should change the names (see the [README](pawsey_slurm/README.md) for a simple command). You might also have single reads that are Illumina sequencing, in which case, you should use the [pawsey minion](pawsey_minion/) pipeline to process that data.

See the verions:
   - [pawsey slurm](pawsey_slurm) -- use this for paired end (R1 and R2) reads. Although it's designed to run on Pawsey's setonix, it will probably work on any system with a `/scratch` drive
   - [pawsey minion](pawsey_minion) -- use this for single end reads. Also designed to run on Pawsey's setonix.
   - [slurm](slurm/README.md) - designed to work on Flinders deepthought infrastructure. This is esoteric and probably not portable.
   - [pbs](pbs/README.md) - designed to work on the NCI infrastructrue. 

Currently the pipeline depends on these software

   - fastp
   - samtools
   - minimap2
   - megahit
   - mmseqs
   - vamb

You can install all of these with:

```
mamba env create -f ~/atavide_lite/atavide_lite.yaml
```

Note: if you are using an ephemeral system like Pawsey, we also have a mechanism for making temporary conda installations. See  the [pawsey slurm](pawsey_slurm/README.md) or [pawsey minion](pawsey_minion/README.md) READMEs.


If you use atavide light, please [cite it](citation.cff) and then please also [cite the other papers](references.bib) that describe these great tools.

