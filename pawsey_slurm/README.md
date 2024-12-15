# Process a metagenome through Rob's 2023 pipeline

This is the same scripts as in [../slurm](../slurm) but I've "optimised" this to run on Pawsey. They are not optimised at all, but they run!


# All commands in one go:

```
mkdir -p slurm_output/host_slurm  slurm_output/megahit_slurm  slurm_output/mmseqs_slurm  slurm_output/vamb_slurm slurm_output/fastp_slurm
find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt

export NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')
SRC=~/atavide_lite/slurm
PAWSEY_SRC=~/atavide_lite/pawsey_slurm

HUMANDLDJOB=$(sbatch --parsable $PAWSEY_SRC/download_human.slurm)
JOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 $SRC/fastp.slurm)
HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB --dependency=afterok:$HUMANDLDJOB $SRC/host_removal.slurm)
FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB $SRC/fastq2fasta.slurm)
MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$FAJOB $SRC/mmseqs_easy_taxonomy.slurm)
sbatch --dependency=afterok:$MMSEQSJOB $SRC/mmseqs_taxonomy.slurm
SSJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB --array=1-$NUM_R1_READS:1 $PAWSEY_SRC/mmseqs_add_subsystems_taxonomy.slurm)
sbatch --dependency=afterok:$SSJOB $PAWSEY_SRC/count_subsystems.slurm
MEGAHITJOB=$(sbatch  --parsable --dependency=afterok:$HOSTJOB --array=1-$NUM_R1_READS:1 $SRC/megahit.slurm)
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB $SRC/vamb_concat.slurm)
VMJOB=$(sbatch --parsable  --dependency=afterok:$VCJOB --array=1-$NUM_R1_READS:1 $SRC/vamb_minimap.slurm)
VAMBJOB=$(sbatch --parsable --dependency=afterany:$VMJOB $SRC/vamb.slurm)
CHECKMJOB=$(sbatch --parsable --dependency=afterany:$VAMBJOB $SRC/checkm.slurm vamb/bins/ vamb/checkm)

```



