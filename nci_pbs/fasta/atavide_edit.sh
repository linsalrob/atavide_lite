
# A bash script to run atavide_lite on a single R1 and R2 file

##################################################################
#                                                                #
# Step 1. Install the software we need!                          #
#                                                                #
##################################################################



mamba install -y fastp samtools minimap2 megahit mmseqs2

# we create a new environment for VAMB because it has some odd python dependencies
mamba create -yn vamb vamb


##################################################################
#                                                                #
# Step 2. Define the R1 and R2 files that we are going to        #
#         process. We will use these variables many times        #
#                                                                #
#                                                                #
#         Please put these files in a directory called fastq     #
#                                                                #
##################################################################


R1=22FHV26_S42_R1_001.fastq.gz
R2=22FHV26_S42_R2_001.fastq.gz

if [[ ! -e fastq/$R1 ]]; then
	echo "ERROR: File fastq/$R1 does not exist. Please check if it is there" >&2
	exit 1;
fi

if [[ ! -e fastq/$R2 ]]; then
	echo "ERROR: File fastq/$R2 does not exist. Please check if it is there" >&2
	exit 1;
fi

echo "We are processing fastq/$R1 and fastq/$R2"

##################################################################
#                                                                #
#  Step 3. Run fastp to trim adapters and remove low quailty     #
#          sequencces.                                           #
#                                                                #
##################################################################

mkdir --parents fastq_fastp
fastp -n 1 -l 100 -i fastq/$R1 -I fastq/$R2 -o fastq_fastp/$R1 -O fastq_fastp/$R2 --adapter_fasta ~/GitHubs/atavide_lite/adapters/IlluminaAdapters.fa --thread 8 

##################################################################
#                                                                #
#  Step 4. Use minimap2 to map the fastq files to the human      #
#          genome. We use the human genome from NCBI             #
#          see: https://linsalrob.github.io/ComputationalGenomicsManual/Deconseq #
#                                                                #
##################################################################

mkdir human_mapped
minimap2 -t 8  --split-prefix=tmp$$ -a -xsr  ~/Downloads/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz fastq_fastp/$R1 fastq_fastp/$R2  | samtools view -bh - > human_mapped/output.bam

##################################################################
#                                                                #
#  Step 5. Separate human and non-human data using samtools      #
#                                                                #
##################################################################

mkdir unmapped_reads
samtools fastq -F 3584 -f 77 human_mapped/output.bam  | gzip -c > unmapped_reads/$R1
samtools fastq -F 3584 -f 141 human_mapped/output.bam  | gzip -c > unmapped_reads/$R2 
samtools fastq -f 4 -F 1 human_mapped/output.bam  | gzip -c > unmapped_reads/singletons.fastq.gz 


##################################################################
#                                                                #
#  Step 6. Assemble the unmapped reads, including singletons     #
#                                                                #
##################################################################

megahit -1 unmapped_reads/$R1 -2 unmapped_reads/$R2 -r unmapped_reads/singletons.fastq.gz -o megahit


##################################################################
#                                                                #
#  Step 6 (alt). Assemble unmapped reads, including singletons   #
#                As an alternative to using megahit, you can     #
#                use spades. We find this sometimes doesn't work #
#                well with lots of metagenomes, but usually      #
#                works well with one or two fastq files like     #
#                this example.                                   #
#                                                                #
##################################################################

spades.py -1 no_human/$R1 -2 no_human/$R2 -o spades_assembly -t 8


##################################################################
#                                                                #
#  Step 7. Use MMSeqs to compare to the UniRef50 database        #
#          You can download the UniRef50 database from UniProt   #
#          but I prefer to use the mmseqs2 automatic download    #
#          which includes the preformatted taxonomy files.       #
#                                                                #
#          Here we use mmseqs easy taxonomy, which still gives   #
#          us other output formats we use too.                   #
#                                                                #
#          But, mmseqs2 requires fasta format for easy taxonomy  #
#          so we convert the files to fasta files first!         #
#                                                                #
##################################################################

# If you need to download the UniRef50 databases, you can do this
mkdir UniRef50                                        #
mmseqs2  mmseqs databases --threads 8 UniRef50 UniRef50/UniRef50 $(mktemp -d -p tmp)

# this is a special bash construct that means in the variable $R1 replace "fastq.gz" with ".fasta"
# so our fastq file (22FHV26_S42_R1_001.fastq.gz) becomes (22FHV26_S42_R1_001.fasta). 
R1FASTA=${R1/fastq.gz/fasta}
R2FASTA=${R2/fastq.gz/fasta}

atavide_lite/bin/fastq2fasta unmapped/$R1 fasta/$R1FASTA
atavide_lite/bin/fastq2fasta unmapped/$R2 fasta/$R2FASTA

# the bash construct $(mktemp -d -p tmp) means make a directory in tmp
mmseqs easy-taxonomy fasta/$R1FASTA fasta/$R2FASTA UniRef50/UniRef50 mmseqs $(mktemp -d -p tmp) --start-sens 1 --sens-steps 3 -s 7 --threads 32

# this gives us a taxonomy report, top hit report, and other useful outputs.


