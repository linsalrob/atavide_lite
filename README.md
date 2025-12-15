[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.sdsu.edu/research)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/atavide_lite)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17934071.svg)](https://doi.org/10.5281/zenodo.17934071)

# atavide lite

Atavide lite is a simple, yet complete workflow for metagenomics data analysis, including QC/QA, optional host 
removal, annotation, assembly and cross-assembly, and individual read based annotations. 

The motivation is based on the more complete [atavide pipeline](https://github.com/linsalrob/atavide) we built that 
uses snakemake as a workflow manager. We found that solution to be effective, but routine failure at different steps
was hard to debug and follow. In addition, as we move between compute resources, we need to adjust the time and 
memory requirements for each step, and that was not easy to do with snakemake.

Our goal is to provide a simple, easy to use, and easy to understand pipeline that can be used for metagenomics data 
analysis, but one that is broken down into individual steps that can be run one at a time, and that can be easily
modified to suit your compute resources. 

Our solution is to craft a series of scripts, suitable for different clusters. We still lean on 
[snakemake](https://snakemake.readthedocs.io/en/stable/) for some parts of the pipeline, but we run each
part separately, so it is straightforward to see what has worked and what has failed. 
We provide generic slurm and pbs scripts
that will run on most clusters, and then there are specific scripts for the
[Pawsey Supercomputing Centre](https://pawsey.org.au/) setonix system,
Flinders University's [deepthought](https://doi.org/10.25957/FLINDERS.HPC.DEEPTHOUGHT) cluster, and
the [National Computational Infrastructure](https://nci.org.au/) Gadi system. These are
the machines that we use everyday, and so we maintain those scripts to ensure that they work for us. If you would like
to ammend the scripts for your cluster, please submit a pull request, or open an issue, and we will try to address it.

In our experience, each cluster has enough [minor differences](#system-nuances) that it is easier
to maintain individual scripts for each cluster, rather than trying to make a single script that works on all clusters.

# Pipeline steps

Our pipeline is designed to be run in a series of steps, and each step can be run independently. In our day to day work
we lean heavily on the `--dependency` option in `sbatch` to ensure that each step is run only after the previous step 
has completed successfully.

1. Run [fastp](https://github.com/OpenGene/fastp) to trim Illumina or Nanopore barcodes. We provide those in [adapters](adapters/)
2. Use [minimap2](https://github.com/lh3/minimap2) and [samtools](https://www.htslib.org/) to filter out host and 
not host reads. Currently, host reads are ignored. Host can be human, sharks, coral, or anything else.
3. Use [mmseqs](https://github.com/soedinglab/MMseqs2) 
[easy-taxonomy](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-assignment) to compare not host reads to 
[UniRef](https://www.uniprot.org/). By default we use [UniRef50](https://www.uniprot.org/help/uniref) but you 
could use any version.
4. Create a taxonomic summary for each sample and make a single `.tsv` file.
5. Connect in the subsystems from [BV-BRC](https://www.bv-brc.org/), and make a table that includes 
subsystems and taxonomic information
6. Create a subsystems taxonomy for the data
7. Use [vamb](https://github.com/RasmussenLab/vamb) to bin the reads into MAGs

We have described all the steps in a [detailed description of the workflows](DETAILED_PROCESSING_STEPS.md)

# different versions

## Paired vs Single End

In our current processing, we have:
 - Paired end reads from MGI or Illumina sequencing. Those files usually end `_R1.fastq.gz` and `_R2.fastq.gz`, 
and the code looks for those.
 - Single end reads from ONT sequencing. These files end `.fastq.gz`

We have two versions of the pipeline that work with _either_ paired end or single end, and you need to choose
the appropriate version for your data.

**However:** If you download some sequences from SRA, ENA, or DDBJ they may have paired end reads that 
end `_1.fastq.gz` and `_2.fastq.gz` in which case you should change the names 
(see the [README](pawsey_shortread/README.md) for a simple command). 
You might also have Illumina sequencing single end reads (which is oldschool!), in which case, 
you should use the [pawsey minion](pawsey_minion/README.md) pipeline to process that data.

See the verions:
   - [pawsey shortread](pawsey_shortread) -- use this for paired end (R1 and R2) reads. Although it's designed to run on 
Pawsey's setonix, it will probably work on any system with a `/scratch` drive
   - [pawsey minion](pawsey_minion) -- use this for single end reads. Also designed to run on Pawsey's setonix.
   - [deepthought shortread](deepthought_shortread/README.md) - designed to work on Flinders deepthought infrastructure. This is esoteric 
and probably not portable, because the deepthought system has a `$BGFS` drive that is used for temporary storage.
   - [nci_pbs](nci_pbs/README.md) - designed to work on the NCI infrastructrue. 

Currently the pipeline depends on these software

   - samtools>=1.20
   - fastp
   - minimap2>=2.29
   - checkm-genome
   - mmseqs2
   - megahit
   - rclone
   - rsync
   - parallel
   - pigz
   - pytaxonkit
   - snakemake
   - sra-tools
   - snakemake-executor-plugin-cluster-generic
   - taxonkit

You can install all of these with:

```
mamba env create -f ~/atavide_lite/atavide_lite.yaml
```

Note: if you are using an ephemeral system like Pawsey, we also have a mechanism for making temporary conda installations. See  the [pawsey shortread](pawsey_shortread/README.md) or [pawsey minion](pawsey_minion/README.md) READMEs.


If you use atavide light, please [cite it](citation.cff) and then please also [cite the other papers](references.bib) that describe these great tools.

<a id='system-nuances'></a>
# System nuances

This is not a comprehensive list of the nuances of each system, but it provides some of the differences that we 
run into and the motivation for some of the choices we made in the scripts.

## Flinders' [deepthought](https://doi.org/10.25957/FLINDERS.HPC.DEEPTHOUGHT)

Deepthought uses slurm for scheduling and has a `$BGFS` drive that is used for temporary storage. This is fast local access
that is only available to the compute node that you are running on. It is not shared between nodes. For most
processing, it is a lot quicker to transfer the data to the `$BGFS` drive, and then run the processing there, and 
finally copying the files back to the working directory when they are complete. This is especially true for
processes that use a lot of memory and create temporary files, such as `megahit` and `mmseqs2`.

## Pawsey Supercomputing Centre's [Setonix](https://pawsey.org.au/)

Setonix uses slurm for scheduling and has a `/scratch` drive that is used for temporary storage, 
but is available from the compute nodes.

## NCI's [Gadi](https://nci.org.au/)

Gadi uses PBS for scheduling a `/g/data` drive that is used for temporary storage. Gadi does not allow array jobs, so 
we have to run each step separately.

Gadi also has fast local drives that are accessible from compute nodes, and they are at `$PBS_JOBFS`
