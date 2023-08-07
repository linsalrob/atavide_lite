[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.sdsu.edu/research)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/atavide_lite)

# atavide lite

A simpler version of atavide that relies only on slurm or PBS scripts. Some of the settings maybe specific for our compute resources

[atavide](https://github.com/linsalrob/atavide) is a full-fledged metagenome processing pipeline. This is a simpler version that uses either slurm scripts or PBS scripts to accomplish the same things. In this initial release, we avoided snakemake to work through the steps one at a time, although that will likely come soon.

See the verions:
   - [slurm](slurm/README.md) - designed to work on Flinders deepthought infrastructure
   - [pbs](pbs/README.md) - designed to work on the NCI infrastructrue for now


Currently the pipeline depends on

   - fastp
   - samtools
   - minimap2
   - megahit
   - mmseqs
   - vamb

You can install all of these with:

```
mamba install -y fastp samtools minimap2 megahit mmseqs vamb
```
