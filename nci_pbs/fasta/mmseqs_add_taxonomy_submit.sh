#!/bin/bash

find $(find mmseqs -mindepth 1 -type d -not -exec sh -c 'ls -1 "{}"|egrep -i -q "_tophit_report_subsystems.gz$"' \; -print ) -type f -name \*_tophit_report.gz -exec qsub -v IN={} atavide_lite/mmseqs_add_taxonomy.pbs \;

