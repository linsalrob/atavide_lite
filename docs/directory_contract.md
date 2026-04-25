# Directory Contract: Pipeline Input/Output Documentation

This document defines the inputs, outputs, and validation criteria for each step of the atavide_lite pipeline. Understanding this "contract" helps with debugging, partial runs, and verifying successful completion.

---

## General Conventions

- **Directory naming**: Most steps create a subdirectory named after the step (e.g., `fastq_fastp/`, `mmseqs/`, `vamb/`)
- **File naming**: Output files typically preserve the input sample identifier with step-specific suffixes
- **Compressed files**: Most intermediate and final files are gzip compressed (`.gz`)
- **Success validation**: Check for output file existence, non-zero size, and step-specific metrics

---

## Step 1: Quality Control with fastp

**Purpose**: Trim adapters and low-quality bases from raw sequencing reads.

### Inputs
- **Directory**: `fastq/` (or `SOURCE` from config)
- **Files**: 
  - Paired-end: `*_R1.fastq.gz`, `*_R2.fastq.gz`
  - Single-end: `*.fastq.gz`
- **Adapter file**: `~/atavide_lite/adapters/IlluminaAdapters.fa` (or appropriate)

### Outputs
- **Directory**: `fastq_fastp/`
- **Trimmed reads**: Same naming as input (e.g., `sample_R1.fastq.gz`, `sample_R2.fastq.gz`)
- **Reports**: 
  - `fastp_output/sample_R1.json` - JSON metrics
  - `fastp_output/sample_R1.html` - HTML report

### Validation
- Output files exist and size > 0
- Check JSON for metrics:
  - `summary.after_filtering.total_reads` > 0
  - `summary.after_filtering.q30_rate` (quality metric)
- Typical reduction: 5-15% of reads filtered

### Typical Resources
- **Memory**: 8-32 GB
- **Threads**: 16
- **Time**: 10 min - 2 hours (depending on file size)

### Failure Symptoms
- **Empty output**: Check adapter file path, input file corruption
- **Excessive filtering**: Review quality thresholds (-n, -l parameters)
- **Logs**: Check `slurm_output/fastp_slurm/fastp-*.err`

---

## Step 2: Host Removal with minimap2/samtools

**Purpose**: Separate host and non-host (microbial) reads.

### Inputs
- **Directory**: `fastq_fastp/`
- **Host reference**: Path specified in `HOSTFILE` variable
- **Files**: Trimmed paired or single-end reads

### Outputs
- **Directories**: 
  - `${HOST}/` - reads mapping to host
  - `${HOSTREMOVED}/` - reads NOT mapping to host (microbial)
- **Files**: `sample_R1.fastq.gz`, `sample_R2.fastq.gz` in each directory
- **BAM files** (optional): `sample.bam`, `sample.bam.bai`

### Validation
- Both directories exist with reads
- Check read counts: `zcat file.fastq.gz | echo $(($(wc -l)/4))`
- Sum of host + non-host should ≈ input reads (allowing for unmapped)
- Non-host directory should contain microbial reads (typically 60-99% of total)

### Typical Resources
- **Memory**: 32-64 GB (depends on reference genome size)
- **Threads**: 16-32
- **Time**: 30 min - 4 hours

### Failure Symptoms
- **No non-host reads**: Wrong host reference, or sample is pure host
- **All reads assigned to host**: Check minimap2 parameters, reference index
- **Logs**: Check stderr for minimap2/samtools errors

---

## Step 3: Convert fastq to fasta

**Purpose**: Prepare sequences for mmseqs2 (requires fasta format).

### Inputs
- **Directory**: `fastq_fastp/` or `${HOSTREMOVED}/`
- **Files**: `*_R1.fastq.gz`, `*_R2.fastq.gz`

### Outputs
- **Directory**: `fasta/`
- **Files**: `sample_R1.fasta.gz`, `sample_R2.fasta.gz`
- Headers should include `/1` or `/2` suffix for paired reads

### Validation
- Fasta files exist and size > 0
- Record count matches fastq: `zgrep -c '^>' file.fasta.gz` should equal fastq reads
- Headers are properly formatted (no spaces, /1 or /2 suffix)

### Typical Resources
- **Memory**: 4-16 GB
- **Threads**: 1
- **Time**: 5-30 minutes

### Failure Symptoms
- **Empty output**: Input file corruption, conversion tool issue
- **Header format errors**: Will cause mmseqs2 to fail in next step

---

## Step 4: Taxonomic Classification with mmseqs2 easy-taxonomy

**Purpose**: Assign taxonomy to reads using protein database (e.g., UniRef50).

### Inputs
- **Directory**: `fasta/`
- **Database**: `${MMSEQS_DB_DIR}/${MMSEQS_DB_NAME}`
- **Files**: Paired fasta files

### Outputs
- **Directory**: `mmseqs/sample/`
- **Files**:
  - `sample_tophit_aln.gz` - alignment results
  - `sample_tophit_report.gz` - taxonomic assignments
  - `sample_lca.tsv.gz` - lowest common ancestor assignments
- **Logs**: mmseqs creates detailed log files

### Validation
- Output directory exists with multiple files
- Check report file size > 0
- Count classified reads: `zcat sample_tophit_report.gz | wc -l`
- Typical classification rate: 30-70% of reads

### Typical Resources
- **Memory**: 128-220 GB (database-dependent)
- **Threads**: 32-64
- **Time**: 2-24 hours (highly variable)

### Failure Symptoms
- **Database errors**: Check database path and integrity
- **Out of memory**: Increase allocation or use smaller database
- **Temp directory full**: Specify larger `TMPDIR`
- **Logs**: Check `slurm_output/mmseqs_slurm/mmseqsET-*.err`

---

## Step 5: Taxonomic Summarization

**Purpose**: Aggregate taxonomic results across samples and create summary tables.

### Inputs
- **Directory**: `mmseqs/sample/`
- **Files**: `sample_tophit_report.gz` from all samples

### Outputs
- **Directory**: `taxonomy_summary/` or similar
- **Files**: 
  - `all_samples_taxonomy.tsv` - combined taxonomy table
  - `sample_counts.tsv` - per-sample read counts by taxon
  - Optional: visualizations (e.g., Sankey plots)

### Validation
- Summary files exist and size > 0
- Row counts match number of unique taxa
- Counts sum correctly across samples

### Typical Resources
- **Memory**: 8-32 GB
- **Threads**: 1-8
- **Time**: 10 min - 2 hours

---

## Step 6: Functional Annotation (BV-BRC Subsystems)

**Purpose**: Connect taxonomic assignments to functional subsystems.

### Inputs
- **Taxonomy results**: `mmseqs/sample/sample_tophit_report.gz`
- **BV-BRC mappings**: Database files in `${BVBRC_DIR}`

### Outputs
- **Directory**: `subsystems/` or within `mmseqs/`
- **Files**: 
  - `sample_subsystems_taxonomy.tsv` - integrated table
  - `subsystems_summary.tsv` - functional pathway counts

### Validation
- Subsystems file exists with taxonomic and functional columns
- Check for expected subsystem categories (e.g., Metabolism, Cell Wall)
- Count annotated reads vs. classified reads

### Typical Resources
- **Memory**: 16-64 GB
- **Threads**: 4-16
- **Time**: 30 min - 4 hours

---

## Step 7: Assembly with MEGAHIT

**Purpose**: Assemble reads into contigs for downstream binning.

### Inputs
- **Directory**: `fastq_fastp/` or `${HOSTREMOVED}/`
- **Files**: Paired or single-end reads for assembly

### Outputs
- **Directory**: `megahit/sample/`
- **Files**:
  - `final.contigs.fa` - assembled contigs
  - `log` - assembly log with statistics
  - `intermediate_contigs/` - k-mer intermediate results

### Validation
- `final.contigs.fa` exists and size > 0
- Check log for assembly statistics:
  - Number of contigs
  - N50 (typically 500-5000 bp for metagenomes)
  - Total assembled length
- Contig count: `grep -c '^>' final.contigs.fa`

### Typical Resources
- **Memory**: 64-256 GB (highly sample-dependent)
- **Threads**: 16-32
- **Time**: 1-48 hours (size and complexity-dependent)

### Failure Symptoms
- **Out of memory**: Reduce k-mer range or use more memory
- **No contigs**: Low coverage, poor quality reads
- **Logs**: Check stderr and `megahit/sample/log`

---

## Step 8: Binning with VAMB

**Purpose**: Bin assembled contigs into metagenome-assembled genomes (MAGs).

### Inputs
- **Contigs**: `megahit/sample/final.contigs.fa` from multiple samples
- **Concatenated contigs**: `vamb_concat/contigs_concat.fna`
- **BAM files**: Alignment of reads to contigs

### Outputs
- **Directory**: `vamb/`
- **Files**:
  - `clusters.tsv` - contig-to-bin assignments
  - `vae_clusters.tsv` - alternative clustering
  - `bins/` - directory with individual MAG fasta files
  - `latent.npz` - VAE latent representation
  - `lengths.npz` - contig lengths
  - `log.txt` - VAMB run log

### Validation
- `clusters.tsv` exists with assignments: `clustername<TAB>contigname`
- Bin count: `cut -f1 clusters.tsv | sort -u | wc -l`
- Check bins directory for individual MAG files
- Validate MAG completeness with CheckM (separate step)

### Typical Resources
- **Memory**: 64-128 GB
- **Threads**: 8-16 (GPU optional)
- **Time**: 2-24 hours

### Failure Symptoms
- **Version incompatibility**: Check VAMB version and parameter names
- **Empty bins**: Low coverage, insufficient samples, or poor assembly
- **CUDA errors**: GPU issues (VAMB works on CPU too)
- **Logs**: Check `vamb/log.txt` and stderr

---

## Quick Validation Script Template

```bash
#!/bin/bash
# Quick validation for atavide_lite pipeline outputs

SAMPLE="your_sample_name"

echo "=== Step 1: fastp ==="
[ -f "fastq_fastp/${SAMPLE}_R1.fastq.gz" ] && echo "✓ Trimmed reads exist" || echo "✗ Missing trimmed reads"
[ -f "fastp_output/${SAMPLE}_R1.json" ] && echo "✓ QC report exists" || echo "✗ Missing QC report"

echo "=== Step 2: Host removal ==="
[ -d "${HOSTREMOVED}" ] && echo "✓ Non-host directory exists" || echo "✗ Missing non-host directory"
[ -s "${HOSTREMOVED}/${SAMPLE}_R1.fastq.gz" ] && echo "✓ Non-host reads exist" || echo "✗ Empty non-host reads"

echo "=== Step 3: Fasta conversion ==="
[ -f "fasta/${SAMPLE}_R1.fasta.gz" ] && echo "✓ Fasta files exist" || echo "✗ Missing fasta files"

echo "=== Step 4: mmseqs2 ==="
[ -d "mmseqs/${SAMPLE}" ] && echo "✓ mmseqs directory exists" || echo "✗ Missing mmseqs directory"
[ -f "mmseqs/${SAMPLE}/${SAMPLE}_tophit_report.gz" ] && echo "✓ Taxonomy report exists" || echo "✗ Missing report"

echo "=== Step 7: Assembly ==="
[ -f "megahit/${SAMPLE}/final.contigs.fa" ] && echo "✓ Assembly exists" || echo "✗ Missing assembly"

echo "=== Step 8: VAMB ==="
[ -f "vamb/clusters.tsv" ] && echo "✓ VAMB clusters exist" || echo "✗ Missing clusters"
```

---

## Notes

- **Partial runs**: You can restart from any step if previous outputs exist
- **Debugging**: Always check both stdout and stderr logs
- **Resource tuning**: These are guidelines; adjust based on your data size and cluster
- **Parallel processing**: Many steps use array jobs; track individual job success
