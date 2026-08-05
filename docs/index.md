# atavide lite documentation

`atavide_lite` is a modular metagenomics workflow for quality control, optional host-read removal, read-based taxonomic and functional profiling, assembly, and recovery of metagenome-assembled genomes (MAGs). It supports paired-end short reads and single-end or long reads on Slurm- and PBS-managed high-performance computing (HPC) systems.

The project deliberately uses visible, independently submitted batch scripts instead of hiding the entire analysis inside one workflow invocation. That makes failures easier to locate, individual stages easier to repeat, and CPU, memory, wall-time, storage, and GPU requests easier to tune for each cluster.

## Where to begin

New users should follow this path:

1. Read [Overview](overview.md) to decide whether the workflow fits the analysis.
2. Follow [Installation](installation.md) to clone the repository, build its utilities, and create the software environments.
3. Use [Choosing a profile](profiles.md) to select scripts for the read type and HPC system.
4. Complete [Quick start](quickstart.md) with a small, non-sensitive dataset.
5. Refer to [Configuration](configuration.md) and [Pipeline stages](pipeline.md) while scaling up.

Cluster maintainers and contributors can start with [Contributing](contributing.md) and [Adding a cluster](cluster-support.md).

```{toctree}
:maxdepth: 2
:caption: Introduction

overview
design
```

```{toctree}
:maxdepth: 2
:caption: User guide

installation
profiles
quickstart
configuration
pipeline
outputs
troubleshooting
```

```{toctree}
:maxdepth: 2
:caption: Technical guide

methods
utilities
cluster-support
contributing
```

## Project links

- [Source repository](https://github.com/linsalrob/atavide_lite)
- [Issue tracker](https://github.com/linsalrob/atavide_lite/issues)
- [License](https://github.com/linsalrob/atavide_lite/blob/main/LICENSE)
- [Citation metadata](https://github.com/linsalrob/atavide_lite/blob/main/citation.cff)

The documentation describes the repository as shipped. Cluster policies and third-party software change independently, so confirm scheduler limits and software versions with the relevant HPC provider before submitting expensive jobs.
