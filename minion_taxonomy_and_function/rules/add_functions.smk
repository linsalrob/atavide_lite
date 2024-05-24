"""

Add the functions to the mmseqs output

"""

import os
import sys


# NOTE: Will update this once I have the data provisioned on ROADs, but my allocation is
# full, so I requested more :)

rule get_sqlite:
    output:
        sql = os.path.join(SCRATCH, "uniref.sqlite")
    shell:
        """
        # curl -Lo {output.sql} http://
        rclone copy A18:dbs/uniref.sqlite {output.sql}
        """

rule add_ss_to_tax:
    input:
        sql = os.path.join(SCRATCH, "uniref.sqlite"),
        thr = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report.gz")
    output:
        ss = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report_subsystems.gz")
    script:
        "../scripts/easy_taxonomy_to_function.py"

