#!/bin/bash

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

: "${PAWSEY_PROJECT:?PAWSEY_PROJECT is not set}"

ENV="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"

if [[ ! -d "$ENV" ]]; then
    echo "ERROR: environment directory does not exist: $ENV" >&2
    exit 1
fi

if [[ ! -x "$ENV/bin/python" ]]; then
    echo "ERROR: $ENV does not look like a valid conda env: missing bin/python" >&2
    exit 1
fi

if [[ ! -d "$ENV/conda-meta" ]]; then
    echo "ERROR: $ENV does not look like a valid conda env: missing conda-meta/" >&2
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available on PATH" >&2
    exit 1
fi

if ! conda run -p "$ENV" --no-capture-output python -c 'import sys; print("Conda prefix:", sys.prefix)'; then
    echo "ERROR: conda cannot run commands inside $ENV" >&2
    exit 1
fi

for exe in samtools fastp minimap2 mmseqs megahit rclone rsync parallel pigz snakemake fasterq-dump taxonkit; do
    if ! conda run -p "$ENV" bash -lc "command -v $exe >/dev/null"; then
        echo "ERROR: missing executable in environment: $exe" >&2
        exit 1
    fi
done

echo "OK: $ENV is our atavide_lite conda environment. Continuing!"
exit 0
