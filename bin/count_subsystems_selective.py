"""
Colate the subsystems output from the mmseqs output into single files for different levels.

In this version, we will remove any of the taxa provided, so we neeed
to at least have the taxa added.

Now the format is:
0       UniRef50_A0A5J4P836
1       12
2       0.954
3       7.245
4       0.653
5       85831
6       species
7       Bacteroides acidifaciens
8       Pyruvate,phosphate dikinase (EC 2.7.9.1)                8 --> Function
9       Energy                                                  9 --> Class
10      Energy and Precursor Metabolites Generation            10 --> Level 1
11      Central Metabolism                                     11 --> Level 2
12      Glycolysis and Gluconeogenesis                         12 --> Subsystem
13      4.0                                                    13 --> Weighted count
14      k__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides acidifaciens --> taxonomy
        

Note that we also check that the taxa column is the weight_column + 1

"""

import gzip
import os
import re
import sys
import argparse
from atavide_lib import colours
from collections import defaultdict

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', '--directory', help='directory of mmseqs outputs [mmseqs]', default='mmseqs')
    parser.add_argument('-s', '--subsystems', help='subsystems output directory [subsystems]', default='subsystems')
    parser.add_argument('-m', '--metadata', help='metadata file. If you provide this we strip of _S34 (or whatever) and then try and rename your columns')
    parser.add_argument('-t', '--taxa', help='taxa to remove. You can supply multiple taxa. e.g. p__Bacteroidota, f__Pseudomonadaceae', action='append', required=False, default=[])
    parser.add_argument('-n', '--name', help='sample name for the file outputs', default='atavide')


    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.directory):
        print(f"FATAL: {args.directory} doesm not exist. Did you specify the right mmseqs output directory?", file=sys.stderr)
        sys.exit(1)

    metadata = {}
    sub = re.compile(r'_S\d+$')
    if args.metadata:
        sys.stderr.write(f"Reading metadata from {args.metadata}\n")
        with open(args.metadata, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                metadata[p[0]] = sub.sub('', p[1])
                metadata[sub.sub('', p[0])] = p[1]

    all_samples = set()
    total = {}
    ss_total = defaultdict(dict)
    ss_class = defaultdict(dict)
    all_classes = set()
    ss_lvl1 = defaultdict(dict)
    all_lvl1 = set()
    ss_lvl2 = defaultdict(dict)
    all_lvl2 = set()
    ss_sub = defaultdict(dict)
    all_sub = set()
    ss_all = defaultdict(dict)
    all_subsystems = set()
    ss_fn = defaultdict(dict)
    all_fn = set()

    weight_column = None
    taxa_column = None
    for sample in os.listdir(args.directory):

        # mmseqs/SAGCFN_22_00789_S34/SAGCFN_22_00789_S34_tophit_report_subsystems.gz

        ss_file = os.path.join(args.directory, sample, f"{sample}_tophit_report_subsystems_taxa.gz")

        if not os.path.exists(ss_file):
            sys.stderr.write(f"Skipping {sample} because there is no subsystems file with the taxa added\n")
            continue

        if args.verbose:
            print(f"Reading {sample} from {ss_file}", file=sys.stderr)
        sample_id =  sub.sub('', sample)
        if sample_id in metadata:
            sample_id = metadata[sample_id]

        all_samples.add(sample_id)
        total[sample_id] = 0
        with gzip.open(ss_file, 'rt') as f:
            for lcount,l in enumerate(f):
                p = l.strip().split("\t")
                if len(p) < 15:
                    print(f"{colours.RED}Error: {sample} has {len(p)} columns, not 15 or 16{colours.ENDC}", file=sys.stderr)
                    sys.exit(-1)
                # check the input file format
                if lcount == 0:
                    for idx in [13, 14]:
                        try:
                            float(p[idx])
                            real_weight_col = idx
                            break
                        except ValueError:
                            pass
                    if weight_column and weight_column != real_weight_col:
                        print(f"{colours.PINK}WARNING: Originally we were using {weight_column} as the weights, and now"
                              f"we are using {real_weight_col} for the weights. Did we mix taxa formats??", file=sys.stderr)
                        weight_column = real_weight_col
                    elif not weight_column:
                        weight_column = real_weight_col
                    # check the taxa column is where we expect it to be
                    taxa_column = weight_column + 1
                    if 'k__' not in p[taxa_column]:
                        print(f"{colours.RED}Error: We expected the taxa column to be column {taxa_column} based on the weight column, but it doesn't look like a taxa column. Please check your input files{colours.ENDC}", file=sys.stderr)
                        sys.exit(-1)
                    if args.verbose:
                        print(f"{colours.GREEN}Using {real_weight_col} as the weight column and {taxa_column} for taxonomy {colours.ENDC}",
                              file=sys.stderr)


                total[sample_id] = total.get(sample_id, 0) + float(p[weight_column]) ## this is the total of all reads
                # skip this read if it is in one of the taxa we want to remove
                toskip = False
                for taxon in args.taxa:
                    if taxon in p[taxa_column]:
                        toskip = True
                        break
                if toskip:
                    continue

                # if we don't have p[9] --> top level, we don't count ss.
                if p[9]:
                    myval = float(p[weight_column]) 
                    ss_total[sample_id] = ss_total.get(sample_id, 0) + myval ## this is the total of only those reads that have a subsystems match

                    if sample_id not in ss_class:
                        ss_class[sample_id] = {}
                    ss_class[sample_id][p[9]] = ss_class[sample_id].get(p[9], 0) + myval  
                    all_classes.add(p[9])

                    if sample_id not in ss_lvl1:
                        ss_lvl1[sample_id] = {}
                    ss_lvl1[sample_id][p[10]] = ss_lvl1[sample_id].get(p[10], 0) + myval
                    all_lvl1.add(p[10])

                    if sample_id not in ss_lvl2:
                        ss_lvl2[sample_id] = {}
                    ss_lvl2[sample_id][f"{p[10]}; {p[11]}"] = ss_lvl2[sample_id].get(f"{p[10]}; {p[11]}", 0) + myval
                    all_lvl2.add(f"{p[10]}; {p[11]}")

                    if sample_id not in ss_sub:
                        ss_sub[sample_id] = {}
                    ss_sub[sample_id][p[12]] = ss_sub[sample_id].get(p[12], 0) + myval
                    all_sub.add(p[12])

                    if sample_id not in ss_fn:
                        ss_fn[sample_id] = {}
                    all_with_fn = f"{p[9]}; {p[10]}; {p[11]}; {p[12]}; {p[8]}"
                    ss_fn[sample_id][all_with_fn] = ss_fn[sample_id].get(all_with_fn, 0) + myval
                    all_fn.add(all_with_fn)

                    if sample_id not in ss_all:
                        ss_all[sample_id] = {}
                    all_with_sub = f"{p[9]}; {p[10]}; {p[11]}; {p[12]}"
                    ss_all[sample_id][all_with_sub] = ss_all[sample_id].get(all_with_sub, 0) + myval
                    all_subsystems.add(all_with_sub)

    # now we have all the data, let's write it out
    os.makedirs(args.subsystems, exist_ok=True)

    ## raw counts
    if args.verbose:
        print('Writing raw counts', file=sys.stderr)

    sorted_samples = sorted(all_samples)
    with open(f"{args.subsystems}/{args.name}_class_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_classes):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_class[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_level1_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl1):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl1[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_level2_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl2):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl2[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_subsystems_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_sub):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_sub[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_all_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_subsystems):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_all[sample].get(ss, 0)))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_all_fn_raw.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_fn):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_fn[sample].get(ss, 0)))
            out.write("\n")



    ## normalized counts

    ## All normalizations
    #
    # This normalisation is to the total number of reads that have any match in the database

    if args.verbose:
        print('Writing normalised counts (all)', file=sys.stderr)

    for sample in total:
        total[sample] /= 1e6

    with open(f"{args.subsystems}/{args.name}_class_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_classes):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_class[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_level1_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl1):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl1[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_level2_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl2):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl2[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_subsystems_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_sub):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_sub[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_all_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_subsystems):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_all[sample].get(ss, 0) / total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_all_fn_norm_all.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_fn):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_fn[sample].get(ss, 0) / total[sample]))
            out.write("\n")



    ## Subsystems normalization

    if args.verbose:
        print('Writing normalised counts (subsystems)', file=sys.stderr)

    # we divide the totals by 1e6 so that when we do the division we now are normalized per million reads
    for sample in ss_total:
        ss_total[sample] /= 1e6

    with open(f"{args.subsystems}/{args.name}_class_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_classes):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_class[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_level1_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl1):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl1[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_level2_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_lvl2):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_lvl2[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_subsystems_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_sub):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_sub[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_all_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_subsystems):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_all[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")

    with open(f"{args.subsystems}/{args.name}_all_fn_norm_ss.tsv", 'w') as out:
        out.write("\t" + "\t".join(sorted_samples) + "\n")
        for ss in sorted(all_fn):
            out.write(ss)
            for sample in sorted_samples:
                out.write("\t" + str(ss_fn[sample].get(ss, 0) / ss_total[sample]))
            out.write("\n")



    with open(f"{args.subsystems}/README.md", 'w') as out:
        out.write(f"""
# NORMALIZATIONS

Currently we perform three normalizations:

1. {args.name}_*_raw.tsv

This is the non-normalised data, so just the raw counts. For each sequence, if it appears in one subsystem we incremenet that count by 1, but if it occurs in more than one subsystem, we increment that count by 1/n (1/2 for 2 subsystems, 1/3 for 3 subsystems, etc).

2. {args.name}_*_norm_all.tsv

This data is normalised for _all_ reads, regardless of whether they are in a subsystem or not. This makes smaller numbers. 

3. {args.name}_*_norm_ss.tsv

This data is normalised only to the number of reads that match to subsystems, so if there is a lot of other stuff we ignore it.

""")

