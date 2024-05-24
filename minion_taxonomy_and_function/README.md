
# MinION Metagenome Taxonomy and Function

This pipeline calculates the read-based taxonomy and function for your metagenomes.

Before you start, you need to run these three steps separately:

1. Dorado – basecalling
2. Filtlong – QC/QA – removes low quality (<5%) and v. short reads (<1,000bp) (You might also check for adapter sequences with fastp, but that should be removed by Dorado)
3. Host removal with minimap2 and samtools filters

This pipeline then predicts _all_ open reading frames, and calculates their taxonomy and function

4. ORF calling – direct on fastq – >30 amino acids. This will generate multiple orfs per sequence
5. mmseqs2 is used to generate easy taxonomy on amino acids
6. We add taxonomy labels using taxonkit
7. Summarise taxonomy per contig (or long sequence)
8. Make one table each for kingdom, phylum, class, order, family, genus, species
9. We add functions and subsystems from the BVBRC
10. We summarise those subsystems (see below).


# Taxonomy

We provide several taxonomy level files:

- kingdom.tsv.gz
- phylum.tsv.gz
- class.tsv.gz
- order.tsv.gz
- family.tsv.gz
- genus.tsv.gz
- species.tsv.gz



# Functions

We add the subsystems (when known) to the UniProt functions, and also normalise the count of reads across proteins with multiple functions (e.g. if protein 1 is in two subsystems each gets 0.5 * read count).

The columns in the output are

Column | Example | Meaning
--- | --- | ---
0 |  UniRef50_L1JF06 | Protein ID
1 |  1 | Number of reads that map
2 |  0.636 | err...
3 |  0.636 | err ...
4 |  0.338 | err...
5 |  131567 | taxonomy ID
6 |  no rank | rank
7 |  cellular organisms | Taxonomy description 
8 |  Phytoene synthase (EC 2.5.1.32) | Protein function
9 |  Metabolism | Subsystem Class
10 | Fatty Acids, Lipids, and Isoprenoids | Subsystem Level 1
11 | Steroids and Hopanoids | Subsystem Level 2
12 | Hopanoid biosynthesis | Subsystem 
13 | Phytoene synthase (EC 2.5.1.32) | Function
14 | 0.5 | Weighted count

