#!/bin/bash
# Task file: sample<TAB>chunk_root<TAB>merge_work_directory<TAB>canonical_output_directory

#SBATCH --job-name=mmseqs_merge
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=slurm_output/mmseqs_slurm/merge-%A_%a.out
#SBATCH --error=slurm_output/mmseqs_slurm/merge-%A_%a.err

set -euo pipefail
TASKS=${MERGE_TASKS:?Set MERGE_TASKS to a tab-delimited task file.}
task_index=${SLURM_ARRAY_TASK_ID:?Submit as an array.}
line=$(sed -n "${task_index}p" "$TASKS")
IFS=$'\t' read -r sample chunk_root merged_output output_dir <<< "$line"
[[ -n "$sample" && -d "$chunk_root" && -n "$merged_output" && -n "$output_dir" ]] || {
    echo 'ERROR: invalid merge task.' >&2
    exit 2
}
marker="$merged_output/MERGE.DONE"
if [[ -s "$marker" ]]; then
    for suffix in lca.tsv.gz report.gz tophit_aln.gz tophit_report.gz; do
        [[ -s "$output_dir/${sample}_${suffix}" ]] || {
            echo "ERROR: $marker exists but ${sample}_${suffix} is missing" >&2
            exit 2
        }
    done
    echo "$sample already merged"
    exit 0
fi

MMSEQS=${MMSEQS_BIN:-/scratch/${PAWSEY_PROJECT}/${USER}/software/miniconda3/atavide_lite/bin/mmseqs}
TARGET=${TARGET_DB:-/scratch/${PAWSEY_PROJECT}/${USER}/Databases/UniRef100/UniRef100}
[[ -x "$MMSEQS" && -s "${TARGET}.dbtype" ]] || {
    echo 'ERROR: MMseqs or target database is unavailable.' >&2
    exit 2
}
SCRIPT_DIR=${MMSEQS_CHUNK_MERGE_DIR:-$HOME/GitHubs/atavide_lite/pawsey_minion/mmseqs_chunk_merge}
[[ -s "$SCRIPT_DIR/restore_chunk_taxonomy.py" && -s "$SCRIPT_DIR/validate_tophit_report.py" ]] || {
    echo "ERROR: MMseqs chunk/merge module unavailable: $SCRIPT_DIR" >&2
    exit 2
}

mapfile -t chunk_dirs < <(find "$chunk_root" -mindepth 1 -maxdepth 1 -type d -name 'chunk_*' | sort)
(( ${#chunk_dirs[@]} >= 1 )) || { echo 'ERROR: no chunk directories found.' >&2; exit 2; }

lca_files=()
aggregate_reports=()
alignments=()
tophit_reports=()
result_dbs=()
swapped_dbs=()
for chunk_dir in "${chunk_dirs[@]}"; do
    [[ -s "$chunk_dir/DONE" ]] || { echo "ERROR: incomplete chunk: $chunk_dir" >&2; exit 2; }
    mapfile -t found_lca < <(find "$chunk_dir" -mindepth 1 -maxdepth 1 -type f -name '*_lca.tsv')
    mapfile -t found_report < <(find "$chunk_dir" -mindepth 1 -maxdepth 1 -type f -name '*_report')
    mapfile -t found_aln < <(find "$chunk_dir" -mindepth 1 -maxdepth 1 -type f -name '*_tophit_aln')
    mapfile -t found_tophit < <(find "$chunk_dir" -mindepth 1 -maxdepth 1 -type f -name '*_tophit_report.gz')
    mapfile -t found_result < <(find "$chunk_dir/tmp" -type f -name result.dbtype -printf '%h/result\n')
    mapfile -t found_swapped < <(find "$chunk_dir/tmp" -type f -name result_aln_swapped.dbtype -printf '%h/result_aln_swapped\n')
    if (( ${#found_lca[@]} != 1 || ${#found_report[@]} != 1 || ${#found_aln[@]} != 1
          || ${#found_tophit[@]} != 1 || ${#found_result[@]} != 1 || ${#found_swapped[@]} != 1 )); then
        echo "ERROR: expected one complete easy-taxonomy output set in $chunk_dir" >&2
        exit 2
    fi
    lca_files+=("${found_lca[0]}")
    aggregate_reports+=("${found_report[0]}")
    alignments+=("${found_aln[0]}")
    tophit_reports+=("${found_tophit[0]}")
    result_dbs+=("${found_result[0]}")
    swapped_dbs+=("${found_swapped[0]}")
done

mkdir -p "$merged_output" "$output_dir"
work="$merged_output/tmp.${SLURM_JOB_ID:-manual}.${task_index}"
[[ ! -e "$work" ]] || { echo "ERROR: temporary path already exists: $work" >&2; exit 2; }
mkdir "$work"
trap 'rm -rf -- "$work"' EXIT
threads=${SLURM_CPUS_PER_TASK:-1}

echo "$(date -Is) $sample: concatenating per-query LCA and top-hit alignment outputs"
lca_count_file="$work/lca_rows"
awk -F '\t' -v count_file="$lca_count_file" '
    NF != 8 { printf "ERROR: LCA line %d has %d fields\n", NR, NF > "/dev/stderr"; bad=1; exit 2 }
    $1 in seen { printf "ERROR: duplicate LCA query: %s\n", $1 > "/dev/stderr"; bad=1; exit 2 }
    { seen[$1]=1; print }
    END { if (!bad) print NR > count_file }
' "${lca_files[@]}" | pigz -p "$threads" > "$work/${sample}_lca.tsv.gz"
lca_rows=$(<"$lca_count_file")

alignment_count_file="$work/alignment_rows"
awk -F '\t' -v count_file="$alignment_count_file" '
    NF != 12 { printf "ERROR: alignment line %d has %d fields\n", NR, NF > "/dev/stderr"; bad=1; exit 2 }
    { print }
    END { if (!bad) print NR > count_file }
' "${alignments[@]}" | pigz -p "$threads" > "$work/${sample}_tophit_aln.gz"
alignment_rows=$(<"$alignment_count_file")

echo "$(date -Is) $sample: regenerating aggregate taxonomy report"
if (( ${#result_dbs[@]} == 1 )); then
    pigz -p "$threads" -c "${aggregate_reports[0]}" > "$work/${sample}_report.gz"
else
    current=${result_dbs[0]}
    for (( i=1; i<${#result_dbs[@]}; i++ )); do
        combined="$work/taxonomy_result_${i}"
        "$MMSEQS" concatdbs "$current" "${result_dbs[$i]}" "$combined"
        current=$combined
    done
    "$MMSEQS" taxonomyreport "$TARGET" "$current" "$work/${sample}_report" \
        --report-mode 0 --threads "$threads"
    pigz -p "$threads" "$work/${sample}_report"
fi
report_top_level=$(pigz -dc "$work/${sample}_report.gz" | awk -F '\t' '
    NF != 6 { printf "ERROR: report line %d has %d fields\n", NR, NF > "/dev/stderr"; exit 2 }
    $5 == 0 || $5 == 1 { total += $2 }
    END { print total + 0 }
')
(( report_top_level == lca_rows )) || {
    echo "ERROR: aggregate report top-level count ($report_top_level) != LCA rows ($lca_rows)" >&2
    exit 2
}

echo "$(date -Is) $sample: regenerating and validating top-hit report"
if (( ${#swapped_dbs[@]} == 1 )); then
    cp "${tophit_reports[0]}" "$work/${sample}_tophit_report.gz"
else
    "$MMSEQS" mergedbs "$TARGET" "$work/result_aln_swapped" "${swapped_dbs[@]}"
    "$MMSEQS" summarizealis "$work/result_aln_swapped" "$work/result_aln_swapped_sum" --threads "$threads"
    "$MMSEQS" addtaxonomy "$TARGET" "$work/result_aln_swapped_sum" "$work/result_aln_swapped_sum_tax" --threads "$threads"
    "$MMSEQS" createtsv "$TARGET" "$work/result_aln_swapped_sum_tax" "$work/${sample}_tophit_report" --threads "$threads"
    pigz -p "$threads" "$work/${sample}_tophit_report"
    restore_args=(--merged-report "$work/${sample}_tophit_report.gz" --output "$work/${sample}_tophit_report.corrected.gz")
    for report in "${tophit_reports[@]}"; do restore_args+=(--chunk-report "$report"); done
    /usr/bin/python3.11 "$SCRIPT_DIR/restore_chunk_taxonomy.py" "${restore_args[@]}"
    mv "$work/${sample}_tophit_report.corrected.gz" "$work/${sample}_tophit_report.gz"
fi
validation_args=(--merged-report "$work/${sample}_tophit_report.gz" --output "$work/tophit_validation.json")
for report in "${tophit_reports[@]}"; do validation_args+=(--chunk-report "$report"); done
/usr/bin/python3.11 "$SCRIPT_DIR/validate_tophit_report.py" "${validation_args[@]}"

for path in "$work/${sample}_lca.tsv.gz" "$work/${sample}_report.gz" \
            "$work/${sample}_tophit_aln.gz" "$work/${sample}_tophit_report.gz"; do
    pigz -t "$path"
done
mv -f "$work/${sample}_lca.tsv.gz" "$output_dir/${sample}_lca.tsv.gz"
mv -f "$work/${sample}_report.gz" "$output_dir/${sample}_report.gz"
mv -f "$work/${sample}_tophit_aln.gz" "$output_dir/${sample}_tophit_aln.gz"
mv -f "$work/${sample}_tophit_report.gz" "$output_dir/${sample}_tophit_report.gz"
mv "$work/tophit_validation.json" "$merged_output/tophit_validation.json"
printf 'sample=%s\nchunks=%s\nlca_rows=%s\ntophit_alignment_rows=%s\nreport_top_level=%s\n' \
    "$sample" "${#chunk_dirs[@]}" "$lca_rows" "$alignment_rows" "$report_top_level" > "$marker"
echo "$(date -Is) $sample: all canonical MMseqs outputs completed and validated"
