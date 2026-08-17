# Process a metagenome through Rob's atavide lite pipeline

`atavide lite` is a series of slurm scripts that run sequentially.  Notice that in the bulk of the script runs, we use the `--dependency=afterok:` to check if the previous job completed correctly. If the previous job does not complete properly, the later jobs do not run. That means that this pipelines completes each step before the next one runs, and it is immmediately apparent when something has crashed. It is not often apparent _why_ something has crashed, but that is a problem for you to solve!

These scripts are adapted to run on Flinder's deepthought computer. They are designed to work with single end long read sequencing

The differences between this version and others are:

1. Since we don't have R1 and R2 files, things are simpler. We only need a `reads.txt` file with the read names.
2. We use slightly different `fastp` parameters.
3. The download scripts will touch the files if they are already present
4. As noted elsewhere, deepthought uses a fast local file system called `$BGFS` that we copy data on and off.

## 00. Create a conda installation or refresh your existing one.

On Pawsey, its best to make the conda installation fresh to ensure all the deleted files are replaced!

```
$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh 2> /dev/null  &&  \
    mamba env remove --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite && \
    mamba env remove --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite_vamb
mamba env create --yes --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite --file ../atavide_lite.yaml
mamba env create --yes --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite_vamb --file ../atavide_lite_vamb.yaml
```

Make sure you build atavide_lite/bin:

```
pushd $HOME/GitHubs/atavide_lite/bin
make all
popd
```

## 0. Get your data.

We expect you to start with a directory called `fastq` that has your reads in it.

## 1. Make a SRC variable that will make life easier for you

This is just a shortcut to make life easier for you!

```
export SRC=$HOME/GitHubs/atavide_lite/pawsey_minion
```

## 1. Create a DEFINITIONS.sh file

Copy _one_ of the DEFINITIONS files to the same location as your `fastq` directory, and renamed it to DEFINITIONS.sh. 

We have several examples:
    - DEFINITIONS_human.sh - removes human genomes
    - DEFINITIONS_shark.sh - removes shark genomes
    - DEFINITIONS_mouse.sh - removes mouse genomes
    - DEFINITIONS_nohost.sh - does not do any host genome removal

If you have a different host to remove, I'm sure you can use one of those as a template!

We also have a `sample name` variable in that file, which we use in some of the processing.

For example, to remove human DNA you can do:

```
cp $SRC/DEFINITIONS_human.sh DEFINITIONS.sh
nano DEFINITIONS.sh
```


# 2. Check the names of your files.

We need to be consistent with the names of our files, so we start by looking at the file names.

Mostly with long read sequencing, we only remove the `.fastq.gz` part of the name, and you leave the rest. 


## 3. Create a reads.txt file and a NUM_READS variable

```
find fastq -type f -printf "%f\n" > reads.txt
export NUM_READS=$(wc -l reads.txt | cut -f 1 -d ' ')
echo "There are $NUM_READS samples to process"
```

## 4. Create the output directorys for the slurm scripts.

This just keeps the output tidy!

```
mkdir -p slurm_output/host_slurm  slurm_output/megahit_slurm  slurm_output/mmseqs_slurm  slurm_output/vamb_slurm slurm_output/fastp_slurm
```

## 5. Download some databases

We are nearly ready to process the data, but we need some databases like the human genome, the UniRef100 database, and so on. We 
can download those simultaneously using the cluster:

```
HUMANDLDJOB=$(sbatch --parsable $PAWSEY_SRC/download_human.slurm)
TAXDLDJOB=$(sbatch --parsable $PAWSEY_SRC/download_taxon_db.slurm)
UNIREFJOB=$(sbatch --parsable $PAWSEY_SRC/download_uniref50.slurm)
UNIREF100JOB=$(sbatch --parsable $PAWSEY_SRC/download_uniref100.slurm)
VAMB_INSTALL=$(sbatch --parsable $PAWSEY_SRC/vamb_install.slurm)
```

I would wait until those databases have downloaded before continuing, however you can also do the first (quality control) step because it doesn't require a database.

If you are using a different host genome, of course omit the human genome download (the first slurm command here), and download your genome instead.



## 6. Quailty control of the data

We use `fastp` for quality control of the data:

```
JOB=$(sbatch --parsable --array=1-$NUM_READS:1 $SRC/fastp.slurm)
```

## 7. Remove the host DNA

`host_removal` uses whatever is defined in the DEFINITIONS.sh file, be it human, mouse, shark, or ...

```
HOSTJOB=$(sbatch --parsable --array=1-$NUM_READS:1 --dependency=afterok:$JOB $SRC/host_removal.slurm)
```

## 8. Assemble the host removed sequences

We use these assemblies, e.g. for binning with VAMB. The assemblies are not used in the rest of the read based annotations, so you don't need to wait for this
to complete before submitting the subsequent jobs.

```
MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --array=1-$NUM_READS:1 $SRC/megahit.slurm)
```

## 9. Convert to fasta.

Please switch to the `atavide_lite` directory and make the executables before you run this step.

```
pushd $HOME/GitHubs/atavide_lite/bin
make all
popd
```

```
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB $SRC/fastq2fasta.slurm)
```

## 10. Run mmseqs.

Note, that we are using the UniRef 100 which gives more hits than UniRef50. Please read [this comparison of UniRef50 vs UniRef100](https://fame.flinders.edu.au/blog/2026/05/24/uniprot)

For UniRef100:

```
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_READS:1 --dependency=afterok:$FAJOB $SRC/mmseqs_easy_taxonomy.slurm)
```

## 11. Summarise the taxonomy from the mmseqs output files

This creates the `taxonomy` and `taxonomy_summary` directories with tables of taxa.

```
sbatch --dependency=afterok:$MMSEQSJOB $SRC/mmseqs_summarise_taxonomy.slurm
```

## 12. Add the subsystems and count them

This creates the `subsystems` directory with tables of subsystem counts

```
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_READS:1 $SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
COUNTSSJOB=$(sbatch --parsable --dependency=afterok:$SSJOB $SRC/count_subsystems.slurm)
```

## 13. Create the data for a sankey plot or other comparisons.

I like the Sankey plots to visualise where the data went, and also to check that the analysis doesn't throw up unexpected errors. Take a look in the `sankey_plot.txt` output file
created by this command for directions on how to make the figure.

```
SANKEYJOB=$(sbatch --parsable --dependency=afterok:$COUNTSSJOB $SRC/sankey_plot.slurm)
```

## VAMB for binning

If you have the assembly from above, you can use VAMB for binning the contigs

```
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_READS:1 $SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB $SRC/vamb.slurm)
CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB $SRC/checkm.slurm vamb/bins/ vamb/checkm)
```



# All commands in one go:

## create a conda environment
```
mamba create --yes atavide_lite --file ../atavide_lite.yaml
mamba activate atavide_lite
export ATAVIDE_CONDA=atavide_lite
```

## Run commands

```
mkdir -p slurm_output/host_slurm  slurm_output/megahit_slurm  slurm_output/mmseqs_slurm  slurm_output/vamb_slurm slurm_output/fastp_slurm
find fastq -type f -printf "%f\n" > reads.txt

export NUM_READS=$(wc -l reads.txt | cut -f 1 -d ' ')
echo $NUM_READS

SRC=~/GitHubs/atavide_lite/pawsey_minion
cp $SRC/DEFINITIONS.sh .

# edit the DEFINITIONS file to change the sample name

JOB=$(sbatch --parsable --array=1-$NUM_READS:1 $SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_READS:1 --dependency=afterok:$JOB $SRC/host_removal.slurm)
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB $SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_READS:1 --dependency=afterok:$FAJOB $SRC/mmseqs_easy_taxonomy.slurm)
sbatch --dependency=afterok:$MMSEQSJOB $SRC/mmseqs_summarise_taxonomy.slurm
sbatch --dependency=afterok:$MMSEQSJOB $SRC/read_fate.slurm
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_READS:1 $SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
COUNTSSJOB=$(sbatch --parsable --dependency=afterok:$SSJOB $SRC/count_subsystems.slurm)
SANKEYJOB=$(sbatch --parsable --dependency=afterok:$COUNTSSJOB $SRC/sankey_plot.slurm)

MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --array=1-$NUM_READS:1 $SRC/megahit.slurm)
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_READS:1 $SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB $SRC/vamb.slurm)
CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB $SRC/checkm.slurm vamb/bins/ vamb/checkm)

```
