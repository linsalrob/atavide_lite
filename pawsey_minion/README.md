# Process a metagenome through Rob's 2023 pipeline

This is the same scripts as in [../minion](../minion) but I've "adapted" this to run on Pawsey. They are not optimised at all, but they run!

Notes:
1. The download scripts will touch the files if they are already present


# All commands in one go:

## create a conda environment
```
TMP=$(for i in {1..12}; do printf "%x" $((RANDOM % 16)); done)
mamba env create --yes --prefix /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP --file ../atavide_lite.yaml
mamba activate /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP
export ATAVIDE_CONDA=/scratch/pawsey1018/edwa0468/software/miniconda3/$TMP
```

## Run commands

```
mkdir -p slurm_output/host_slurm  slurm_output/megahit_slurm  slurm_output/mmseqs_slurm  slurm_output/vamb_slurm slurm_output/fastp_slurm
find fastq -type f -printf "%f\n" > reads.txt

export NUM_R1_READS=$(wc -l reads.txt | cut -f 1 -d ' ')
SRC=~/atavide_lite/pawsey_minion
cp $SRC/DEFINITIONS.sh .

# edit the DEFINITIONS file to change the sample name

# download the databases
HUMANDLDJOB=$(sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/download_human.slurm)
MMSEQSDLD=$(sbatach --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_download.slurm)

# note: If you don't need those two, then submit the next line and use MMSEQSDLD=$JOB; HUMANDLDJOB=$JOB

JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB,$HUMANDLDJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/host_removal.slurm)
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB,$MMSEQSDLD --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_easy_taxonomy.slurm)
sbatch --dependency=afterok:$MMSEQSJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_taxonomy.slurm
sbatch --dependency=afterok:$MMSEQSJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/read_fate.slurm
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
sbatch --dependency=afterok:$SSJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/count_subsystems.slurm
MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA --array=1-$NUM_R1_READS:1 $SRC/megahit.slurm)
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb.slurm)
CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/checkm.slurm vamb/bins/ vamb/checkm)

```


These are the jobs without any dependencies, in case you run them one at a time:

```
JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/host_removal.slurm)
FAJOB=$(sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_easy_taxonomy.slurm)
sbatch --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_taxonomy.slurm
sbatch --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/read_fate.slurm
SSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/mmseqs_add_subsystems_taxonomy_fast.slurm)
sbatch --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/count_subsystems.slurm
MEGAHITJOB=$(sbatch  --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA --array=1-$NUM_R1_READS:1 $SRC/megahit.slurm)
VCJOB=$(sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --array=1-$NUM_R1_READS:1 --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/vamb.slurm)
CHECKMJOB=$(sbatch --parsable --export=ATAVIDE_CONDA=$ATAVIDE_CONDA $SRC/checkm.slurm vamb/bins/ vamb/checkm)

```


