




# Functions

We add the subsystems (when known) to the UniProt functions, and also normalise the count of reads across proteins with multiple functions (e.g. if protein 1 is in two subsystems each gets 0.5 * read count).

The columns in the output are

Column | Example | Meaning
--- | --- | ---
0 |  UniRef50_L1JF06 | Protein ID
1 |  1 | Number of reads that map
2 |  0.636 | err...
3 |  0.636 | err ...
4 |  0.338 | err...
5 |  131567 | taxonomy ID
6 |  no rank | rank
7 |  cellular organisms | Taxonomy description 
8 |  Phytoene synthase (EC 2.5.1.32) | Protein function
9 |  Metabolism | Subsystem Class
10 | Fatty Acids, Lipids, and Isoprenoids | Subsystem Level 1
11 | Steroids and Hopanoids | Subsystem Level 2
12 | Hopanoid biosynthesis | Subsystem 
13 | Phytoene synthase (EC 2.5.1.32) | Function
14 | 0.5 | Weighted count

