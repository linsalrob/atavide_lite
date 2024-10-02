# Generate a better taxonomy information for atavide

We're going to use [taxonkit](https://bioinf.shenwei.me/taxonkit/) since they have put the work in to make it work

# taxonkit

## data

You need to download the NCBI data. You should then set the $TAXONKIT_DB environment 

taxonkit only requires names.dmp, nodes.dmp, delnodes.dmp and merged.dmp

## set an environment variable

```
DATE=`date +%Y%m`
TAXONKIT_DB=~/Databases/NCBI/taxonomy/taxonomy.$DATE
mkdir -p $TAXONKIT_DB
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $TAXONKIT_DB
```

# snakemake

You need snakemake version 7 to run this. Don't use version 8 at the moment, because it breaks a lot of things.

If you have a directory `mmseqs`, and in that directory you have a bunch of samples, each one is a directory, and in those you have the output from `mmseqs easy-taxonomy` like this:

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
snakemake --profile slurm -s ~/GitHubs/atavide_lite/new_taxonomy/taxonomy.smk
```

and it should make the `*.lca.taxonomy.tsv.gz` files for you
