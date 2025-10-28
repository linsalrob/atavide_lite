
# A bash script to run atavide_lite on a single R1 and R2 file


##################################################################
#                                                                #
# Define your conda environment and the different                #
# files that you want to use for host, etc here.                 #
#                                                                #
#                                                                #
#                                                                #
##################################################################

ATAVIDE_CONDA=atavide_lite
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

# where is the human genome?
HOSTFILE=$HOME/Databases/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

# where is the mmseqs database for UniRef50? This should be a file and not a directory
UNIREF=$HOME/Databases/UniRef50/UniRef50
# If you need to download the UniRef50 databases, you can do this
mkdir UniRef50
mmseqs2  mmseqs databases --threads 8 UniRef50 UniRef50/UniRef50 $(mktemp -d -p tmp)

# what should we prepend to the output files?
SAMPLENAME=atavide_lite_example

# where is the UniRef50 functions database? This connects UniRef IDs to subsystems
UNIREFFUNC=$HOME/Databases/uniref.sqlite 




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
#  This comes from fastp.slurm                                   #
#                                                                #
##################################################################

mkdir --parents fastq_fastp
fastp -n 1 -l 100 -i fastq/$R1 -I fastq/$R2 -o fastq_fastp/$R1 -O fastq_fastp/$R2 --adapter_fasta $HOME/atavide_lite/adapters/IlluminaAdapters.fa --thread 16

##################################################################
#                                                                #
#  Step 4. Use minimap2 to map the fastq files to the human      #
#          genome. We use the human genome from NCBI             #
#          see: https://linsalrob.github.io/ComputationalGenomicsManual/Deconseq #
#                                                                #
#  This comes from host_removal.slurm                            #
#                                                                #
##################################################################

mkdir human_mapped
minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $HOSTFILE fastq_fastp/$R1 fastq_fastp/$R2 | samtools view -bh | samtools sort -o  human_mapped/output.bam -
samtools fastq -F 3588 -f 65 host_bamfiles/output.bam | gzip -c > human_mapped/$R1
samtools fastq -F 3588 -f 129 host_bamfiles/output.bam | gzip -c > human_mapped/$R2

##################################################################
#                                                                #
#  Step 5. Separate human and non-human data using samtools      #
#                                                                #
#  This comes from host_removal.slurm                            #
#                                                                #
##################################################################

mkdir unmapped_reads
samtools fastq -F 3584 -f 77 human_mapped/output.bam  | gzip -c > unmapped_reads/$R1
samtools fastq -F 3584 -f 141 human_mapped/output.bam  | gzip -c > unmapped_reads/$R2 
samtools fastq -f 4 -F 1 human_mapped/output.bam  | gzip -c > unmapped_reads/singletons.fastq.gz 

##################################################################
#                                                                #
#  Step 6. Convert the unmapped reads to fasta format.           #
#                                                                #
#  This comes from fastq2fasta.slurm                             #
#                                                                #
##################################################################

~/atavide_lite/bin/fastq2fasta -n 1 $R1 - | gzip -c > fasta/${R1/fastq/fasta}
~/atavide_lite/bin/fastq2fasta -n 2 $R2 - | gzip -c > fasta/${R2/fastq/fasta}

##################################################################
#                                                                #
#  Step 7. Use mmseqs2 to map the reads to UniRef50              #
#                                                                #
#  This comes from mmseqs_easy_taxonomy.slurm                    #
#                                                                #
##################################################################

mkdir mmseqs
mmseqs easy-taxonomy fasta/$R1 fasta/$R2 $UNIREF mmseqs/atavide $(mktemp -d -p tmp) --threads 64
# compress the ouput
find mmseqs -type f  | parallel -j 32 gzip

##################################################################
#                                                                #
#  Step 8. Summarise the taxonomy from the mmseqs output         #
#                                                                #
#  This comes from mmseqs_summarise_taxonomy.slurm               #
#                                                                #
##################################################################

snakemake --profile slurm -s ~/GitHubs/atavide_lite/summarise_taxonomy/taxonomy.smk
python ~/GitHubs/atavide_lite/summarise_taxonomy/scripts/join_taxonomies.py -t taxonomy -o taxonomy_summary/ -n $SAMPLENAME

##################################################################
#                                                                #
#  Step 9. add the subsystems omy from the mmseqs output         #
#                                                                #
#  This comes from mmseqs_summarise_taxonomy.slurm               #
#                                                                #
##################################################################

python ~/atavide_lite/bin/easy_taxonomy_to_function_taxa_fast.py -f mmseqs/atavide_tophit_report.gz -d $UNIREFUNC -o mmseqs/atavide_tophit_report_subsystems_taxa
pigz mmseqs/atavide_tophit_report_subsystems_taxa

##################################################################
#                                                                #
#  Step 10. Count the subsystems and create a single file        #
#                                                                #
#  This comes from count_subsystems.slurm                        #
#                                                                #
##################################################################

perl ~/atavide_lite/bin/count_subsystems.pl -d mmseqs -n $SAMPLENAME
find subsystems -type f -exec pigz {} \;

##################################################################
#                                                                #
#  Step 11. Create a beautiful Sankey plot to summarise the      #
#           data visually.                                       #
#                                                                #
#  This comes from sankey_plot.slurm                             #
#                                                                #
##################################################################

python ~/GitHubs/atavide_lite/bin/sankey_plot.py -r R1_reads.txt  --paired -v



##################################################################
#                                                                #
#  Step 12. Assemble the unmapped reads, including singletons     #
#                                                                #
##################################################################

megahit -1 unmapped_reads/$R1 -2 unmapped_reads/$R2 -r unmapped_reads/singletons.fastq.gz -o megahit


##################################################################
#                                                                #
#  Step 12 (alt). Assemble unmapped reads, including singletons  #
#                As an alternative to using megahit, you can     #
#                use spades. We find this sometimes doesn't work #
#                well with lots of metagenomes, but usually      #
#                works well with one or two fastq files like     #
#                this example.                                   #
#                                                                #
##################################################################

spades.py -1 no_human/$R1 -2 no_human/$R2 -o spades_assembly -t 8

