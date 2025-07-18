# Process a metagenome through Rob's 2023 pipeline

This is the same scripts as in [../slurm](../slurm) but I've "optimised" this to run on Pawsey. They are not optimised at all, but they run!


# If you are processing a lot of SRA runs, you will run into 16S libraries. You can screen for them with this command that works on the fastq directory

```
export NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
mkdir -p slurm_output/sixteen_s
sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA  $PAWSEY_SRC/16S_detection_single.slurm
```
grep 'primary mapped' slurm_output/sixteen_s/*out | perl -ne 'm/(\d+\.\d+)\%/; print "$1\n"' | sort -g | (sed -u 10q ; echo ; tail)
grep 'primary mapped' slurm_output/sixteen_s/*out | perl -ne 'm/(\d+\.\d+)\%/; print "$1\n"' | awk '{s+=$1} END {print s/NR}'
```

# Create a new conda environment:

```
TMP=$(for i in {1..12}; do printf "%x" $((RANDOM % 16)); done)
mamba env create --yes --prefix /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP --file ../atavide_lite.yaml
mamba activate /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP
export ATAVIDE_CONDA=/scratch/pawsey1018/edwa0468/software/miniconda3/$TMP
```

# Renaming files

If you download files that have only `_1.fastq.gz` or `_2.fastq.gz` and need to rename them:

```
for F in *_1.fastq.gz; do mv $F ${F/_1/_R1}; done
for F in *_2.fastq.gz; do mv $F ${F/_2/_R2}; done
```


# All commands in one go:

```
mkdir -p slurm_output/host_slurm  slurm_output/megahit_slurm  slurm_output/mmseqs_slurm  slurm_output/vamb_slurm slurm_output/fastp_slurm
find fastq -name \*_R1\* -printf "%f\n" > R1_reads.txt
export NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
echo $NUM_R1_READS

SRC=~/atavide_lite/slurm
PAWSEY_SRC=~/atavide_lite/pawsey_slurm

cp $SRC/DEFINITIONS.sh .

# edit the DEFINITIONS file to change the sample name

HUMANDLDJOB=$(sbatch --parsable $PAWSEY_SRC/download_human.slurm)
TAXDLDJOB=$(sbatch --parsable $PAWSEY_SRC/download_taxon_db.slurm)
UNIREFJOB=$(sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA download_uniref50.slurm)

JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $PAWSEY_SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB --dependency=afterok:$HUMANDLDJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $PAWSEY_SRC/host_removal.slurm)
sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA  $PAWSEY_SRC/16S_detection_single.slurm
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA  $PAWSEY_SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB --dependency=afterok:$UNIREFJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA  $PAWSEY_SRC/mmseqs_easy_taxonomy.slurm)
sbatch --dependency=afterok:$MMSEQSJOB --dependency=afterok:$TAXDLDJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA  $PAWSEY_SRC/mmseqs_taxonomy.slurm
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA   $PAWSEY_SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
COUNTSSJOB=$(sbatch --parsable --dependency=afterok:$SSJOB $PAWSEY_SRC/count_subsystems.slurm)
SANKEYJOB=$(sbatch --parsable --dependency=afterok:$COUNTSSJOB $PAWSEY_SRC/sankey_plot.slurm)

MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --array=1-$NUM_R1_READS:1 $SRC/megahit.slurm)
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_R1_READS:1 $SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB $SRC/vamb.slurm)
CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB $SRC/checkm.slurm vamb/bins/ vamb/checkm)

```


Note: Even with the _fast_ implementation of adding functions to taxonomy, there is still an issue with the code timing out. You can rerun the code
and it will read the last results and then continue from there.

Check for failed jobs:

```
grep CANCELLED slurm_output/mmseqs_slurm/*err
```

Rerun just the failed jobs (in this example, it was #3 and 18 that failed to finish):

```
SSJOB=$(sbatch --parsable  --array=3,18 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA   $PAWSEY_SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
```
