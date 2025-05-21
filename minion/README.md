# Process a metagenome through Rob's 2023 pipeline

These scripts are adapted to run on Flinder's deepthought computer. They are designed to work with single end long read sequencing

Notes:
1. The download scripts will touch the files if they are already present


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

export NUM_R1_READS=$(wc -l reads.txt | cut -f 1 -d ' ')
echo $NUM_R1_READS

SRC=~/atavide_lite/pawsey_minion
cp $SRC/DEFINITIONS.sh .

# edit the DEFINITIONS file to change the sample name

JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 $SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB $SRC/host_removal.slurm)
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB,$MMSEQSDLD --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_easy_taxonomy.slurm)
sbatch --dependency=afterok:$MMSEQSJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_taxonomy.slurm
sbatch --dependency=afterok:$MMSEQSJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/read_fate.slurm
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
COUNTSSJOB=$(sbatch --parsable --dependency=afterok:$SSJOB $SRC/count_subsystems.slurm)
SANKEYJOB=$(sbatch --parsable --dependency=afterok:$COUNTSSJOB $SRC/sankey_plot.slurm)

MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA --array=1-$NUM_R1_READS:1 $SRC/megahit.slurm)
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb.slurm)
CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/checkm.slurm vamb/bins/ vamb/checkm)

```

