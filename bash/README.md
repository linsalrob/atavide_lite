# Atavide Light in a single bash script

Here we present atavide\_light as a single bash script which has all the commands you need.

You probably should not run this as a single script, but run each command one at a time, and make sure that you get the output you expect, and that you understand that output before proceeding onto the next command.

The steps are:

   - Use `fastp` to remove adapters and trim off bad sequences
   - Use `minimap2` to map the reads to the human genome
   - Use `samtools` to separate out the host and non-host reads
   - Use `megahit` or `spades` to assemble the unmapped reads.
   - Use `mmseqs2` to map the reads to UniRef50 to see what is present.

Take a look at the [steps in atavide_light.sh](atavide_light.sh) to see what we do!
