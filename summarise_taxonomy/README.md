# Generate a better taxonomy information for atavide

We're going to use [taxonkit](https://bioinf.shenwei.me/taxonkit/) since they have put the hard work in to make it work well!

# taxonkit

## data

You need to download the NCBI taxonomy data. You should then set the $TAXONKIT_DB environment variable to
point to the source of that data.

`taxonkit` only requires the following files: `names.dmp`, `nodes.dmp`, `delnodes.dmp`, and `merged.dmp`

_Note:_ There was a hiccup when NCBI updated their taxonomy format, but pytaxonkit has a fix for that, 
so you either need to use the latest version of pytaxonkit or to use an older version of the NCBI taxonomy data. 
You can find more details at this [pytaxonkit issue](https://github.com/bioforensics/pytaxonkit/issues/41)

## Set an environment variable for taxonkit

```
DATE=`date +%Y%m`
TAXONKIT_DB=~/Databases/NCBI/taxonomy/taxonomy.$DATE
mkdir -p $TAXONKIT_DB
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $TAXONKIT_DB
```

# Snakemake

This code will run with snakemake version 8 or higher, but you may need to install the 
`snakemake-executor-plugin-cluster-generic` plugin to run it on a cluster (that is included
in the `atavide_lite.yaml` conda environment).

If you have a directory `mmseqs`, and in that directory you have a bunch of samples, each 
one is a directory, and in those you have the output from `mmseqs easy-taxonomy` like this:

```
mmseqs/
    SAGCFN_22_01149_S3/
        SAGCFN_22_01149_S3_lca.tsv.gz    
    SAGCFN_22_01175_S20/
        SAGCFN_22_01175_S20_lca.tsv.gz  
    SAGCFN_23_00214_S1/
        SAGCFN_23_00214_S1_lca.tsv.gz
```

Then you can run the command:

```
snakemake --profile slurm -s ~/GitHubs/atavide_lite/summarise_taxonomy/taxonomy.smk
```

and it will make the `*.lca.taxonomy.tsv.gz` files for you


Once the snakemake has run, you can use this command:

```
python ~/GitHubs/atavide_lite/summarise_taxonomy/scripts/join_taxonomies.py \
   -t taxonomy -o taxonomy_summary/
```

to create the summary file. 

The slurm scripts run both these commands, so you don't need to run them separately.
