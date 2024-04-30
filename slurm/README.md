# Process a metagenome through Rob's 2023 pipeline

This pipeline consists of several separate slurm/perl/python scripts, and we probably want to clean it up at some point. 

For this version, you will need a file `DEFINITIONS.sh` that defines some variables that the scripts use. Notably, the `$FILEEND` variable is used to remove the `.fastq.gz` and other parts of the filename so that you end up with meaningful names. We use this for the R1 filename only.

For example, you might want one of these:

```
FILEEND=_S34_R1.fastq.gz
FILEEND=_S34_L001_R1.fastq.gz
FILEEND=_R1.fastq.gz
```

Currently, we read these variables from `DEFINITIONS.sh`:

```
export FILEEND=_L001_R1_001.fastq.gz
export FAFILEEND=_L001_R1_001.fasta.gz
export SOURCE=fastq
export HOSTREMOVED=fastq_fastp
```

0. Set the source of the scripts

Clone the atavide lite repository to your home directory, and make the compiled code

```
git clone https://github.com/linsalrob/atavide_lite.git
cd atavide_lite/bin
make
SRC=~/atavide_lite/slurm
```

0b. Set up the conda environment

```
mamba env create -f ~/atavide_lite/atavide_lite.yaml
mamba env create -f ~/atavide_lite/atavide_lite_vamb.yaml
```

1. Create a list of R1 files.

Several of the steps will use the file `R1_reads.txt` to know which reads we have. Of course, `snakemake` doesn't need this, but at the moment I'm being very conservative and running each command separately via `slurm`. In addition, this allows me to use `$BGFS` which I don't know how to use with `snakemake`!

```
find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt
NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
echo There are $NUM_R1_READS R1 readsA
if [[ $(find fastq -name \*R2\* | awk 's++{} END {print s}') != $NUM_R1_READS ]]; then echo "There are a different number of R1 and R2 reads"; fi
```


2. Quality control of the sequences


```
JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 $SRC/fastp.slurm)
```

3. HOST REMOVAL (optional)

If you don't want to do host removal just set `HOSTREMOVED=fastq_fastp` in DEFINITIONS.sh

```
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB $SRC/host_removal.slurm)
```

3. Convert to fasta for mmseqs

```
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB $SRC/fastq2fasta.slurm)
```

4. Run mmseqs taxonomy

```
MMSEQSJOB=$(sbatch --parsable --dependency=afterok:$FAJOB $SRC/mmseqs_easy_taxonomy_submit.slurm)
```

4a. Create a taxonomy table

```
sbatch --dependency=afterok:$MMSEQSJOB $SRC/mmseqs_taxonomy.slurm
```

4b. Add the subsystems to the taxonomy

```
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 $SRC/mmseqs_add_subsystems.slurm)
sbatch --dependency=afterok:$SSJOB $SRC/count_subsystems.slurm
```

5. Run megahit
 
```
MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --array=1-$NUM_R1_READS:1 $SRC/megahit.slurm
```

6. Combine assemblies for VAMB

Note: this takes a few (<10) minutes, and so I run it on the short queue

```
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB $SRC/vamb_concat.slurm)
```

7. Map vamb reads.

```
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_R1_READS:1 $SRC/vamb_minimap.slurm)
```


8. Run VAMB

```
sbatch --dependency=afterok:$VMJOB $SRC/vamb.slurm
```


# Grouped VAMB

There comes a time when we want to run vamb on a subset of reads, and so we have a `grouped` version of vamb.
We use a simple `.tsv` file to designate the groups, and then we run each of the components of vamb on that group.

For example, if we have a `.tsv` file with this data:

--- | ---
Host | Sample
Angelshark | GSV317\_
Angelshark | GSV326\_
Angelshark | GSV329\_
Angelshark | GSV342\_
Eagle\_Ray | ER128RecapASuck\_
Eagle\_Ray | ER151Suck\_
Eagle\_Ray | ER247Suck\_


We will separate GSV317, GSV326, GSV329, GSV342 into AngelShark and ER128RecapASuck, ER151Suck, ER247Suck into Eagle\_Ray directories, and then vamb them separately.

The main trick here is to make it so that the sample names only give 0 or 1 results when grepping through R1_reads.txt which is why I append the `_` to these names.

### Step 1: Create the contigs with the groups:

```
sbatch $SRC/vamb_concat_groups.slurm groups_sample.tsv
```

### Step 2: Map all the reads:

```
sbatch --array=1-$NUM_R1_READS:1 $SRC/vamb_minimap_group.slurm
```

### Step 3: Run vamb:

```
sbatch $SRC/vamb_group.slurm
```

All three lines as one go:

```
VCONCAT=$(sbatch --parsable $SRC/vamb_concat_groups.slurm groups_sample.tsv)
VMAP=$(sbatch --parsable  --dependency=afterok:$VCONCAT --array=1-$NUM_R1_READS:1 $SRC/vamb_minimap_group.slurm)
sbatch --parsable  --dependency=afterok:$VMAP $SRC/vamb_group.slurm
```

# All commands in one go:

```
find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt

export NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
SRC=~/atavide_lite/slurm

JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 $SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB $SRC/host_removal.slurm)
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB $SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB $SRC/mmseqs_easy_taxonomy.slurm)
sbatch --dependency=afterok:$MMSEQSJOB $SRC/mmseqs_taxonomy.slurm
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 $SRC/mmseqs_add_subsystems.slurm)
sbatch --dependency=afterok:$SSJOB $SRC/count_subsystems.slurm
MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --array=1-$NUM_R1_READS:1 $SRC/megahit.slurm)
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_R1_READS:1 $SRC/vamb_minimap.slurm)
sbatch --dependency=afterok:$VMJOB $SRC/vamb.slurm


```



