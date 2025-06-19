---
title: 'atavide-lite: simplified atavistic metagenome processing'
tags:
  - Python
  - bash
  - metagenomics
  - dna sequencing
authors:
  - name: Robert A. Edwards
    orcid: 0000-0001-8383-8949
    equal-contrib: false
    affiliation: 1
affiliations:
 - name: Flinders Accelerator for Microbiome Exploration, Flinders University, Bedford Park, South Australia, 5042, Australia
   index: 1
date: 1 July 2025
bibliography: paper.bib
---

# Summary
`atavide lite` is a metagenomics processing pipeline that separates the component parts into
individual steps, allowing users to run only the steps they need. Through experience, we have 
eschewed the use of all-in-one pipelines, but have created a set of modular steps that can be
applied to complex, multi-sample, multi-technology datasets. Together, the atavide-lite processes
combine both read-based and assembly-based approaches to metagenome processing, and can be used
with paired-end or long-read sequencing data. 

# Statement of need

There are many metagenomics pipelines available, and include both web-based pipelines like MG-RAST 
`[@Meyer2008-mw:2008]`, and many, many command-line pipelines 
`[e.g. @Clarke2019-go:2019; @Laudadio2019-og:2019; @Lu2022-yw:2022; @Garfias-Gallegos2022-cy:2022; @Walker2022-yo:2022; 
@Blanco-Miguez2023-ej:2023; @Tyagi2024-wl:2024; @Roach2024-qn:2024]`. Workflow management systems like Nextflow 
`[@Di_Tommaso2017-ir:2017]` and Snakemake `[@Koster2012-qn:2012]` streamlined the creation of pipelines
for bioinformatics analysis, including metagenomics. However, the two major issues with these pipelines
is that they are often fragile, susceptible to breaking when processing thousands of datasets, and 
they struggle to adapt to the nuances of different computational environments.

For example, when processing thousands of samples from the Sequence Read Archive (SRA) 
`[@Leinonen2011-yd:2011]`, we found that some samples failed to download, and we had to 
interrupt the pipeline to fix the issue. Similarly, when using large computations, such as
comparing sequence reads to a database of reference sequences using MMSeqs `[@Steinegger2017-qw:2017]`, 
the computation occasionally times out because of limitations of the compute environment, 
or the computation fails for a variety of other reasons.

Finally, each of the computational platforms we use has nuances for the most
efficient computational processing. For example, our local HPC has a very fast local disk which
is available from a non-standard location and is only available to compute nodes. The Australian 
National Computational Infrastructure's Gadi system does not support large array jobs. The 
Pawsey Supercomputing Centre's Setonix system regularly deletes files on their "/scratch" 
disk, but has a large S3 storage system that is available to all compute nodes.

Therefore, we created a set of modular steps that can be run independently, and that are optimised
for each of the computational platforms we use. We include generic scripts that should be suitable
for any platform, but we encourage users to adapt the scripts to the nuances of their own
computational environment.



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

I acknowledge the support of  everyone

# References
