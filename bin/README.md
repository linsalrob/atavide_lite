# Executables and scripts

These are additional scripts that we need to use for atavide lite. 

Please make the binary files with the following commands before you start using atavide lite. 

```bash
make all
```

Most of the standalone scripts layer additional information into the output files from other code.


# Dependencies

These are some of the dependencies that `atavide\_light` depends upon. You can compile the `C` code by using the [Makefile](Makefile) by typing just `make`.

[vamb_concatenate.py](vamb_concatenate.py) comes from the [VAMB](https://github.com/RasmussenLab/vamb) repository and is provided here for ease of use. Please refer to  [VAMB](https://github.com/RasmussenLab/vamb) and it's documentation, and don't forget to cite [VAMB](https://www.nature.com/articles/s41587-020-00777-4)

Our version also writes `sample_mappings.tsv` next to the output. VAMB renames every contig `S<n>C<contig>`, where `<n>` is the position of its assembly on the command line, and this file records which sample each `S<n>` is:

```
vamb_sample	sample	path
S1	sample_1	/path/to/megahit/sample_1/final.contigs.fa
S2	sample_2	/path/to/megahit/sample_2/final.contigs.fa
```

The sample name comes from the assembly's directory (`megahit/<sample>/final.contigs.fa` or `megahit/<sample>/output/final.contigs.fa`), or from the file name when that is not a generic name such as `final.contigs.fa`. For VAMB output made before this file existed, recreate it with the same inputs, in the same order, and the same `-m`:

```
python vamb_concatenate.py --mapping-only -m 2000 vamb/contigs.fna.gz megahit/*/final.contigs.fa
```

`--mapping-only` checks that each input's first kept contig appears under the expected `S<n>` name in the existing file, and refuses to write the mapping if the order does not match.



