#!/bin/bash
#SBATCH --job-name=mmTax
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --partition=short
#SBATCH -o mmtax-%j.out
#SBATCH -e mmtax-%j.err


set -euo pipefail

python ~/atavide_lite/bin/mmseqs_report_to_table.py -d mmseqs -o mmseqs_taxonomy -v
