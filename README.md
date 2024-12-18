[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.sdsu.edu/research)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/atavide_lite)
[![DOI](https://www.zenodo.org/badge/671497428.svg)](https://www.zenodo.org/badge/latestdoi/671497428)

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

See the verions:
   - [slurm](slurm/README.md) - designed to work on Flinders deepthought infrastructure
   - [pbs](pbs/README.md) - designed to work on the NCI infrastructrue for now
   - [minion pawsey](minion_pawsey)

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


If you use atavide light, please [cite it](citation.cff) and then please also [cite the other papers](references.bib) that describe these great tools.

