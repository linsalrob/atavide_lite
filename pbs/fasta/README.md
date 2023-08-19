# Atavide lite using PBS and fasta files

### Why?

 - we need to use PBS on the [NCI](https://www.nci.org.au/) cluster, hence these scripts. Also, NCI has some ... _oddities_ ... For example, they do not like array jobs, but we can submit lots of individual jobs, and so we work within their preferences
 - we are using fasta because the data was downloaded from the SRA, and so downloading the fasta files is 1/2 the data size. And nobody uses the quality scores anyway!


# Step 0. Separate interleaved fasta files

When you download from SRA using `fasterq-dump` the [quickest option is to ignore the order of the fasta file](https://edwards.flinders.edu.au/fastq-dump/). We separate the fasta files into R1 and R2 reads using [fasta_split](https://github.com/linsalrob/EdwardsLab/blob/master/bin/fasta_split.c).

# Step 1. QC and trimming

We use [cutadapt](https://cutadapt.readthedocs.io/en/stable/) since it is able to handle fasta files. However, it also has limitations: you can't use paired end reads with fasta files

The code [cutadapt_submit.sh](cutadapt_submit.sh) will submit everything in the directory you specify


```
bash ~/GitHubs/atavide_lite/pbs/fasta/cutadapt_submit.sh fasta_split/
```


# Step 2. Remove human contamination

We have a (pbs script](human.pbs) that trims out the R1 and R2 reads, and creates human and non-human directories. Again, there is a submit script to check for R1 and R2 or just R1:


```
bash ~/GitHubs/atavide_lite/pbs/fasta/human_submit.sh cutadapt
```

