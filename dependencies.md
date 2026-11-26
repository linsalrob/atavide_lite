# Software Dependencies

checkm-genome
fastp 
megahit
megahit_toolkit 
minimap2>=2.29
mmseqs2
parallel
pigz
pytaxonkit
rclone
rsync
samtools>=1.20
snakemake
snakemake-executor-plugin-cluster-generic
sra-tools (for analysing SRA data)
taxonkit
vamb 

# Compiled code

There are two compiled applications, 
`~/GitHubs/atavide_lite/bin/fastq2fasta` and 
`~/GitHubs/atavide_lite/bin/fastg2gfa`

You can compile that code with the makefile we've provided.

```
cd ~/GitHubs/atavide_lite/bin/
make
```

