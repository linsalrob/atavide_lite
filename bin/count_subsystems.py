"""
Colate the subsystems output from the mmseqs output into single files for different levels.

Columns in the subsystems file:
 0       UniRef50_L1JF06
 1       1
 2       0.636
 3       0.636
 4       0.338
 5       131567
 6       no rank
 7       cellular organisms
 8       Phytoene synthase (EC 2.5.1.32) --> original function
 9       Metabolism  --> CLASS
 10      Fatty Acids, Lipids, and Isoprenoids --> Level 1
 11      Steroids and Hopanoids --> Level 2
 12      Hopanoid biosynthesis --> Subsystem
 13      Phytoene synthase (EC 2.5.1.32) --> Function
 14      0.5 --> Weighted count

"""
import gzip
import os
import re
import sys
import argparse

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', '--directory', help='directory of mmseqs outputs [mmseqs]', default='mmseqs')
    parser.add_argument('-s', '--subsystems', help='subsystems output directory [subsystems]', default='subsystems')
    parser.add_argument('-m', '--metadata', help='metadata file. If you provide this we strip of _S34 (or whatever) and then try and rename your columns')

    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    metadata = {}
    sub = re.compile(r'_S\d+$')
    if args.metadata:
        sys.stderr.write(f"Reading metadata from {args.metadata}\n")
        with open(args.metadata, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                metadata[p[0]] = sub.sub('', p[1])

    all_samples = set()
    total = {}; ss_total = {}
    ss_class = {}; all_classes = set()
    ss_lvl1 = {}; all_lvl1 = set()
    ss_lvl2 = {}; all_lvl2 = set()
    ss_sub = {}; all_sub = set()
    ss_all = {}; all_subsystems = set()

    for sample in os.listdir(args.directory):

        # mmseqs/SAGCFN_22_00789_S34/SAGCFN_22_00789_S34_tophit_report_subsystems.gz

        if not os.path.exists(os.path.join(args.directory, sample, f"{sample}_tophit_report_subsystems.gz")):
            sys.stderr.write(f"Skipping {sample} because there is no subsystems file\n")
            continue

        if args.verbose:
            print(f"Reading {sample}", file=sys.stderr)
        sample_id =  sub.sub('', sample)
        if sample_id in metadata:
            sample_id = metadata[sample_id]

        all_samples.add(sample_id)
        with gzip.open(os.path.join(args.directory, sample, f"{sample}_tophit_report_subsystems.gz"), 'rt') as f:
            for l in f:
                p = l.strip().split("\t")
                if len(p) != 15:
                    sys.stderr.write(f"Error: {sample} has {len(p)} columns, not 15\n")
                    sys.exit(-1)

                total[sample_id] = total.get(sample_id, 0) + float(p[14]) ## this is the total of all reads
                # if we don't have p[9] --> top level, we don't count ss.
                if p[9]:
                    ss_total[sample_id] = ss_total.get(sample_id, 0) + float(p[14]) ## this is the total of only those reads that have a subsystems match

                    if sample_id not in ss_class:
                        ss_class[sample_id] = {}
                    ss_class[sample_id][p[9]] = ss_class[sample_id].get(p[9], 0) + float(p[14])
                    all_classes.add(p[9])

                    if sample_id not in ss_lvl1:
                        ss_lvl1[sample_id] = {}
                    ss_lvl1[sample_id][p[10]] = ss_lvl1[sample_id].get(p[10], 0) + float(p[14])
                    all_lvl1.add(p[10])

                    if sample_id not in ss_lvl2:
                        ss_lvl2[sample_id] = {}
                    ss_lvl2[sample_id][f"{p[10]}; {p[11]}"] = ss_lvl2[sample_id].get(f"{p[10]}; {p[11]}", 0) + float(p[14])
                    all_lvl2.add(f"{p[10]}; {p[11]}")

                    if sample_id not in ss_sub:
                        ss_sub[sample_id] = {}
                    ss_sub[sample_id][p[12]] = ss_sub[sample_id].get(p[12], 0) + float(p[14])
                    all_sub.add(p[12])

                    if sample_id not in ss_all:
                        ss_all[sample_id] = {}
                    ss_all[sample_id][f"{p[9]}; {p[10]}; {p[11]}; {p[12]}"] = ss_all[sample_id].get(f"{p[9]}; {p[10]}; {p[11]}; {p[12]}", 0) + float(p[14])
                    all_subsystems.add(f"{p[9]}; {p[10]}; {p[11]}; {p[12]}")

    # now we have all the data, let's write it out
    os.makedirs(args.subsystems, exist_ok=True)

    ## raw counts
    if args.verbose:
        print('Writing raw counts', file=sys.stderr)

    sorted_samples = sorted(all_samples)
    with open(f"{args.subsystems}/class_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_classes):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_class[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/level1_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl1):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl1[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/level2_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl2):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl2[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/subsystem_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_sub):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_sub[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/all_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_subsystems):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_all[sample].get(ss, 0)))
            out.write("\n")


    ## normalized counts

    ## All normalizations
    #
    # This normalisation is to the total number of reads that have any match in the database

    if args.verbose:
        print('Writing normalised counts (all)', file=sys.stderr)

    for sample in total:
        total[sample] /= 1e6

    with open(f"{args.subsystems}/class_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_classes):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_class[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/level1_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl1):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl1[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/level2_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl2):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl2[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/subsystem_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_sub):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_sub[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/all_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_subsystems):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_all[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    ## Subsystems normalization

    if args.verbose:
        print('Writing normalised counts (subsystems)', file=sys.stderr)

    # we divide the totals by 1e6 so that when we do the division we now are normalized per million reads
    for sample in ss_total:
        ss_total[sample] /= 1e6

    with open(f"{args.subsystems}/class_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_classes):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_class[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/level1_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl1):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl1[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/level2_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl2):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl2[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/subsystem_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_sub):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_sub[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/all_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_subsystems):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_all[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/README.md", 'w') as out:
        out.write("""
# NORMALIZATIONS

Currently we perform three normalizations:

1. \*\_raw.tsv

This is the non-normalised data, so just the raw counts. For each sequence, if it appears in one subsystem we incremenet that count by 1, but if it occurs in more than one subsystem, we increment that count by 1/n (1/2 for 2 subsystems, 1/3 for 3 subsystems, etc).

2. \*\_norm_all.tsv

This data is normalised for _all_ reads, regardless of whether they are in a subsystem or not. This makes smaller numbers. 

3. \*\_norm_ss.tsv

This data is normalised only to the number of reads that match to subsystems, so if there is a lot of other stuff we ignore it.

""")

