#!/bin/bash
#SBATCH --job-name=write_methods
#SBATCH --time=00:02:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1M
#SBATCH -o %x.out
#SBATCH -e %x.err



if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi


FP_VER=$(fastp -v 2>&1 | cut -d ' ' -f 2)
MM_VER=$(minimap2 --version)
ST_VER=$(samtools --help 2>&1 | grep Version)
TK_VER=$(taxonkit 2>&1 | grep Version)
VB_VER=$(vamb --version 2>/dev/null) || VB_VER="Vamb 5.0.2"

cat <<EOF

Each of the metagenomics sequencing runs was processed separately. In the following commands the variable \$READ is used to indicate a specific run name for long read sequences and \$R1 and \$R2 are used to indicate the R1/R2 paired-end runs.

The long-read data were cleaned with fastp (version $FP_VER) [PMID: 30423086] using the parameters:

fastp -n 2 -l 50 -i fastq/\$READ -o fastq_fastp/\$READ --adapter_fasta ~/atavide_lite/adapters/ONT_RPB114.24.fna

This removes sequences with 2 or more "N" bases, removes the ONT adapters, and requires a minimum read length of 50 bp. This removes low quality Nanopore reads.

The short-read data were cleaned with fastp (version $FP_VER) [PMID: 30423086] using the parameters:

fastp -n 1 -l 100 -i fastq/$R1 -I fastq/$R2 -o fastq_fastp/$R1 -O fastq_fastp/$R2 -j fastp_output/$R1.json -h fastp_output/$R1.html --adapter_fasta atavide_lite/adapters/IlluminaAdapters.fa --thread 16

For the short reads, we remove sequences with a single "N" base, that are shorter than 100 bp, and trim off the Illumina adapters. Those sequences are available in the _atavide_lite_ GitHub repository.

Next, the sequences were compared to the human genome reference sequence GRCh38 using minimap2 (version $MM_VER) [PMID: 29750242] and samtools ($ST_VER) [PMID: 19505943] with the following command for long read sequences:

minimap2 -t 64 --secondary=no --split-prefix=tmp\$\$ -a -x map-ont GRCh38.fna fastq_fastp/\$READ | samtools view -bh | samtools sort -o \$READ.bam -

and the similar command for short read sequences:

minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $HOSTFILE $QC/$R1 $QC/$R2 | samtools view -bh | samtools sort -o host_bamfiles/$R1 -

These commands use the default minimap2 parameters for mapping reads against a reference genome and has been optimised by the tools' authors [PMID: 29750242].

The sequences were identified from the BAM file created by the previous command and filtered to two separate files using samtools  ($ST_VER) [PMID: 19505943]. Sequences that do not match the flag 3588 -- and thus represent mapped reads which pass the platform/vendor quality check, are not PCR or optical duplicates, and are primary alignments -- were filtered as mapping to the human genome. For short reads, additional samtools flags were used to identify R1 reads (matching flag 65) and R2 reads (matching flag 129). Reads that do not match the flag 3584 -- and thus represent unmapped reads which pass the platform/vendor quality check, are not PCR or optical duplicates, and are primary alignments -- are filtered as unmapped reads [PMID: 19505943], and again, short reads were separated using flags that match 77 for R1 reads and 141 for R2 reads.

Reads that did not map to the human genome were mapped against the UniRef50 [PMID: 17379688]  (version database using MMseqs2 (version $MM_VER) with the following command for long reads:

mmseqs easy-taxonomy fasta/\$READS UniRef50 mmseqs/\${READS}_UniRef /tmp --threads 32

and a similar command for short-read sequencing

mmseqs easy-taxonomy fasta/\$R1 fasta/\$R2 UniRef50 mmseqs/\${R1}_UniRef /tmp --threads 32

MMseqs2 easy-taxonomy parameters were optimised by the authors for a balance of sensitivity and speed and include a 7-mer double-match prefilter with compositional bias correction, low-complexity masking with TANTAN, and default similarity thresholds (--k-score ~95) that generate sufficient similar k-mers for accurate detection while maintaining computational efficiency. Candidate hits are further filtered by fast ungapped alignment and refined with vectorized Smith-Waterman alignments, after which taxonomic labels are assigned using a lowest common ancestor (LCA) approach described in more detail in their paper and on their website [PMID: 29035372].

The taxonomies were reformatted with taxonkit [PMID: 34001434] ($TK_VER) and merged into a single table.

The functional annotations were enriched by mapping the UniProt IDs from the MMseq2 output directly to proteins associated with BV-BRC subsystems [PMID: 36350631] (this was a 1:1 mapping and did not require thresholds/cutoffs). The number of reads that MMseqs2 reported mapped to each protein was counted, and the total was divided by the number of mapped reads to normalise the read counts.

For the assembly-based approaches, the metagenomes were assembled using megahit [PMID: 25609793] ($MH_VER) which by default employs a succinct de Bruijn graph approach with multiple k-mer sizes (21 to 141 in steps of 20), automatic memory optimization, and conservative heuristics for contig pruning to balance assembly sensitivity, contiguity, and computational efficiency. Metagenome-assembled genomes (MAGs) were reconstructed using VAMB [PMID: 33398153] ($VB_VER), a variational autoencoder-based approach that jointly encodes k-mer composition and differential coverage profiles of contigs into a low-dimensional latent space. By learning these compressed representations, VAMB effectively separates contigs originating from different microbial genomes and clusters them using a medoid-based algorithm. This deep learning approach has been shown to improve bin purity and completeness compared to conventional binning methods such as MetaBAT [PMID: 31388474] or CONCOCT [PMID: 25218180], particularly in complex metagenomic datasets [PMID: 29183281].


EOF
