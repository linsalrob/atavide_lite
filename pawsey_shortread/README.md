# Process a metagenome through Rob's _atavide lite_ pipeline


`atavide lite` is a series of slurm scripts that run sequentially.  Notice that in the bulk of the script runs, we use the `--dependency=afterok:` to check if the previous 
job completed correctly. If the previous job does not complete properly, the later jobs do not run. That means that this pipelines completes each step before the next one
runs, and it is immmediately apparent when something has crashed. It is not often apparent _why_ something has crashed, but that is a problem for you to solve!

Before you being, please remove any old environments, and create the `atavide_lite` mamba environment, using the command:

```
$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh 2> /dev/null  &&  \
    mamba env remove --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite && \
    mamba env remove --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite_vamb
mamba env create --yes --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite --file ../atavide_lite.yaml
mamba env create --yes --prefix /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite_vamb --file ../atavide_lite_vamb.yaml
```

This creates an `atavide_lite` and `atavide_lite_vamb` conda environments for us to use.


## 0. Get your data.

We expect you to start with a directory called `fastq` that has your reads in it.

## 1. Make a PAWSEY_SRC variable that will make life easier for you

This is just a shortcut to make life easier for you!

```
export PAWSEY_SRC=$HOME/GitHubs/atavide_lite/pawsey_shortread
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
cp $PAWSEY_SRC/DEFINITIONS_human.sh DEFINITIONS.sh
nano DEFINITIONS.sh
```


# 2. Check the names of your files.

We need to be consistent with the names of our files, so we start by looking at the file names.

Sometimes, the files are called `XXXX_R1.fastq.gz` and `XXXX_R2.fastq.gz`, sometimes they are called `XXXX_R1_001.fastq.gz` and `XXXX_R2_001.fastq.gz`, and
sometimes they are called things like `XXXX_L001_R1.fastq.gz` and `XXXX_L0001.R2.fastq.gz`.

Take a look at your fastq files, and decide where you want to trim the name off. Use `nano` to edit the DEFINITIONS.sh file and change the file ending. We will
trim that off all files.

If you download files that have only `_1.fastq.gz` or `_2.fastq.gz` and need to rename them, because we rely on having `R1` and `R2` in the file names.

```
for F in *_1.fastq.gz; do mv $F ${F/_1/_R1}; done
for F in *_2.fastq.gz; do mv $F ${F/_2/_R2}; done
```

## 3. Create an R1_reads.txt file and a NUM_R1_READS variable

Since you have paired end files, we only use the `R1` reads, but use that to get the name of the `R2` reads. The NUM_R1_READS is
used the array jobs where we process each read separately.

```
find fastq -name \*_R1\* -printf "%f\n" > R1_reads.txt
export NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
echo "There are $NUM_R1_READS samples to process"
```

## 4. Create the output directorys for the slurm scripts.

This just keeps the output tidy!

```
mkdir -p slurm_output/host_slurm  slurm_output/megahit_slurm  slurm_output/mmseqs_slurm  slurm_output/vamb_slurm slurm_output/fastp_slurm
```

## 5. Download some databases

We are nearly ready to process the data, but we need some databases like the human genome, the UniRef50 or UniRef100 databases, and so on. We 
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
JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 $PAWSEY_SRC/fastp.slurm)
```

## 7. Remove the host DNA

`host_removal` uses whatever is defined in the DEFINITIONS.sh file, be it human, mouse, shark, or ...

```
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB $PAWSEY_SRC/host_removal.slurm)
```

## 8. Assemble the host removed sequences

We use these assemblies, e.g. for binning with VAMB. The assemblies are not used in the rest of the read based annotations, so you don't need to wait for this
to complete before submitting the subsequent jobs.

```
MEGAHITHR=$(sbatch --parsable --dependency=afterok:$HOSTJOB $PAWSEY_SRC/megahit_hostremoved.slurm)
```

## 9. Convert to fasta.

Please switch to the `atavide_lite` directory and make the executables before you run this step.

```
pushd $HOME/GitHubs/atavide_lite/bin
make all
popd
```

```
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB $PAWSEY_SRC/fastq2fasta.slurm)
```


## 10. Run mmseqs.

Note, that initially I was using UniRef50, however I currently use UniRef 100 which gives more hits. Please read [this comparison of UniRef50 vs UniRef100](https://fame.flinders.edu.au/blog/2026/05/24/uniprot) for more discussion and comparisons.

For UniRef50:

```
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB $PAWSEY_SRC/mmseqs_easy_taxonomy.slurm)
```

For UniRef100:

```
MMSEQS100JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB $PAWSEY_SRC/mmseqs_easy_taxonomy100.slurm)
```

> Note:
> Currently the UniRef100 output is in a directory called `mmseqs_uniref100/mmseqs` so that you can run both UniRef50 and UniRef100. If you only run the UniRef100 analysis 
> (which you probably should) then just move the `mmseqs` directory up to your working directory. e.g. `mv `mmseqs_uniref100/mmseqs ./`

## 11. Summarise the taxonomy from the mmseqs output files

This creates the `taxonomy` and `taxonomy_summary` directories with tables of taxa.

```
sbatch --dependency=afterok:$MMSEQSJOB $PAWSEY_SRC/mmseqs_summarise_taxonomy.slurm
```

## 12. Add the subsystems and count them

This creates the `subsystems` directory with tables of subsystem counts

```
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 $PAWSEY_SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
COUNTSSJOB=$(sbatch --parsable --dependency=afterok:$SSJOB $PAWSEY_SRC/count_subsystems.slurm)
```

## 13. Create the data for a sankey plot or other comparisons.

I like the Sankey plots to visualise where the data went, and also to check that the analysis doesn't throw up unexpected errors. Take a look in the `sankey_plot.txt` output file 
created by this command for directions on how to make the figure.

```
SANKEYJOB=$(sbatch --parsable --dependency=afterok:$COUNTSSJOB $PAWSEY_SRC/sankey_plot.slurm)
```

## 14. Counting 16S genes

If you are processing data from the SRA, you will run undoubtedly into 16S libraries. You can screen for them with this command that works on the fastq directory

```
sbatch --parsable $PAWSEY_SRC/16S_detection_single.slurm
```

Once the counting has finished you can see how many 16S hits there are with:

```
grep 'primary mapped' slurm_output/sixteen_s/*out | perl -ne 'm/(\d+\.\d+)\%/; print "$1\n"' | sort -g | (sed -u 10q ; echo ; tail)
grep 'primary mapped' slurm_output/sixteen_s/*out | perl -ne 'm/(\d+\.\d+)\%/; print "$1\n"' | awk '{s+=$1} END {print s/NR}'
```

## 15. VAMB

Please note: this section of the README needs some improvement, because I have not used the different pieces in a while.

### Using VAMB to bin reads.

Note: if you used the megahit_host_removed the output will be in megahit_all_reads. We cheat, and move that into a directory called `megahit` so we have
`megahit/megahit_all_reads/contigs.fna` as our input file before we continue.

The normal way to use vamb is these three commands:

```
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB $PAWSEY_SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_R1_READS:1 $PAWSEY_SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB --account=${PAWSEY_PROJECT}-gpu $HOME/atavide_lite/pawsey_shortread/vamb.slurm)

CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB $PAWSEY_SRC/checkm.slurm vamb/bins/ vamb/checkm)
```


### Use this code for CROSSASSEMBLY data

<summary>
About the cross-assembly approach
<details>
In this approach, we do a mega-huge assembly with everyting together, and then map the individual reads back and use the complete
contigs with mapped reads in VAMB - it appears that VAMB ignores contigs with no mapped reads, but breaks if you have bam files that 
don't map to any contigs, so we find those and move those out the way.

Then we use the usual VAMB approach to bin the contigs.

</details>
</summary>


Assemble using the `megahit_hostremoved.slurm` script above. This will take a while to run, so do it early! Also note that `megahit` can continue if it is interuppted. Make sure the `--continue` flag is active 
in the `megahit_hostremoved.slurm` script.

```
VCRJOB=$(sbatch --parsable $PAWSEY_SRC/vamb_concat_crass.slurm samples.tsv)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCRJOB --array=1-$NUM_R1_READS:1 $PAWSEY_SRC/vamb_minimap_crass.slurm samples.tsv)
```

#### Sanity check

All the files should be the same length in each directory (but not between directories)

```
find vamb_crass/ -name \*.bam | while read -r BAM; do L=$(samtools view -H $BAM | wc -l); echo -e "$L\t$BAM"; done
```

Then run the rest of vamb:

```
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB --account=${PAWSEY_PROJECT}-gpu $PAWSEY_SRC/vamb_crass.slurm)
sbatch  --dependency=afterok:$VAMBJOB  ~/GitHubs/atavide_lite/pawsey_shortread/vamb_mags_group.slurm samples.tsv
CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB $PAWSEY_SRC/checkm.slurm vamb/bins/ vamb/checkm)
```


## All the commands in one go.

You can, in fact, copy and paste all the commands in one go, once you have the databases installed. If something breaks, it will stop at that point, and you'll need to look in the `slurm` output
directory to find the error.


```
JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 $PAWSEY_SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB $PAWSEY_SRC/host_removal.slurm)
MEGAHITHR=$(sbatch --parsable --dependency=afterok:$HOSTJOB $PAWSEY_SRC/megahit_hostremoved.slurm)
sbatch --parsable $PAWSEY_SRC/16S_detection_single.slurm
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB $PAWSEY_SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB $PAWSEY_SRC/mmseqs_easy_taxonomy.slurm)
MMSEQS100JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB $PAWSEY_SRC/mmseqs_easy_taxonomy100.slurm)
sbatch --dependency=afterok:$MMSEQSJOB $PAWSEY_SRC/mmseqs_summarise_taxonomy.slurm
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 $PAWSEY_SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
COUNTSSJOB=$(sbatch --parsable --dependency=afterok:$SSJOB $PAWSEY_SRC/count_subsystems.slurm)
SANKEYJOB=$(sbatch --parsable --dependency=afterok:$COUNTSSJOB $PAWSEY_SRC/sankey_plot.slurm)
```

## Failed jobs

One of the most common failures is the script taking too long. You can find those, for example if the mmseqs step failed, like this:


```
grep CANCELLED slurm_output/mmseqs_slurm/*err
```

Rerun just the failed jobs (in this example, it was #3 and 18 that failed to finish):

```
SSJOB=$(sbatch --parsable  --array=3,18 $PAWSEY_SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
```
