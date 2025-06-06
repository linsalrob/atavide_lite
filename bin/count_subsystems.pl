use strict;
use Getopt::Std;
use Data::Dumper;

my %opts;
getopts('d:m:s:n:v', \%opts);
unless ($opts{d}) {
	die <<EOF;
	$0
	-d directory of mmseqs outputs [mmseqs]
	-m metadata file. If you provide this we strip of _S34 (or whatever) and then try and rename your columns
	-s subsystems output directory [subsystems]
	-n sample name used for the file outputs
	-v verbose output
EOF
}

unless ($opts{s}) {$opts{s}="subsystems"}
if ($opts{n}) {$opts{n} =~ s/\W//g}
else {$opts{n}="atavide"}


## Input format:

# The original format we used was this, but note that the function is in columns 8 & 13. 
# 0       UniRef50_L1JF06
# 1       1
# 2       0.636
# 3       0.636
# 4       0.338
# 5       131567
# 6       no rank
# 7       cellular organisms
# 8       Phytoene synthase (EC 2.5.1.32)
# 9       Metabolism  --> CLASS
# 10      Fatty Acids, Lipids, and Isoprenoids --> Level 1
# 11      Steroids and Hopanoids --> Level 2
# 12      Hopanoid biosynthesis --> Subsystem 
# 13      Phytoene synthase (EC 2.5.1.32) --> Function
# 14      0.5 --> Weighted count
#
# Then we updated the format and added the taxa in column 15. Then we updated the format again when we added mmseqs_add_subsystems_taxonomy_fast.slurm and
# removed the redundant duplicate of the format.
#
# Now the format is:
# 0       UniRef50_A0A5J4P836
# 1       12
# 2       0.954
# 3       7.245
# 4       0.653
# 5       85831
# 6       species
# 7       Bacteroides acidifaciens
# 8       Pyruvate,phosphate dikinase (EC 2.7.9.1)                8 --> Function
# 9       Energy                                                  9 --> Class
# 10      Energy and Precursor Metabolites Generation            10 --> Level 1
# 11      Central Metabolism                                     11 --> Level 2
# 12      Glycolysis and Gluconeogenesis                         12 --> Subsystem
# 13      4.0                                                    13 --> Weighted count
# 14      k__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides acidifaciens --> taxonomy
#
## This version _only_ works on the 'newer' format from mmseqs_add_subsystems_taxonomy_fast.slurm which has one
# less column (I removed the redundant function column). 
# Its a bit of a mess because it means we have incompatible formats floating around
# However, you can change the format essentially by changing the $weight_column variable

my $weight_column = 13;

# this is a REALLY memory inefficient way of doing this!
my $class; my $ac;
my $lvl1;  my $al1;
my $lvl2;  my $al2;
my $ss;    my $ass;
my $all;   my $aall;

my %allsub;
my %total;
my %sstotal;

my %meta;

if ($opts{m}) {
	open(IN, $opts{m}) || die "$! $opts{m}";
	while (<IN>) {
		chomp;
		my @a=split /\t/;
		$a[0] =~ s/_S\d+$//;
		$meta{$a[0]} = $a[1];
	}
	close IN;
}

opendir(DIR, $opts{d}) || die "$! : $opts{d}";

foreach my $sub (grep {$_ !~ /^\./} readdir(DIR)) {
	# mmseqs/SAGCFN_22_00789_S34/SAGCFN_22_00789_S34_tophit_report_subsystems.gz
	
	my $tax_inputfile = "$opts{d}/$sub/${sub}_tophit_report_subsystems_taxa.gz";
	my $ori_inputfile = "$opts{d}/$sub/${sub}_tophit_report_subsystems.gz";
	my $inputfile;

	if (-e $tax_inputfile) {
		# this is the new format with the taxonomy
		$inputfile = $tax_inputfile;
	} elsif (-e $ori_inputfile) {
		# this is the old format without the taxonomy
		$inputfile = $ori_inputfile;
	} else {
		print STDERR "No subsystems found for $sub (we checked both $tax_inputfile and $ori_inputfile)\n";
		next;
	}

	open(IN, "gunzip -c $inputfile |") || die "$! opening pipe to $inputfile";
	my $id = $sub;
	$id =~ s/_S\d+$//;
	if ($meta{$id}) {$id = $meta{$id}}
	$allsub{$id}=1;
	while (<IN>) {
		chomp;
		my @a=split /\t/;
		$total{$id}+=$a[$weight_column]; ## this is the total of all reads
		next unless ($a[9]);
		$sstotal{$id}+=$a[$weight_column]; ## this is the total of only those reads that have a subsystems match
		$class->{$id}->{$a[9]} += $a[$weight_column];
		$ac->{$a[9]} = 1;
		$lvl1->{$id}->{$a[10]} += $a[$weight_column];
		$al1->{$a[10]} = 1;
		$lvl2->{$id}->{"$a[10]; $a[11]"} += $a[$weight_column];
		$al2->{"$a[10]; $a[11]"} = 1;
		$ss->{$id}->{$a[12]} += $a[$weight_column];
		$ass->{$a[12]} = 1;
		$all->{$id}->{"$a[9]; $a[10]; $a[11]; $a[12]"} += $a[$weight_column];
		$aall->{"$a[9]; $a[10]; $a[11]; $a[12]"} = 1;
	}
	close IN;
}


if ($opts{v}) {
	print STDERR Dumper($all);
}

my @allsub = sort {$a cmp $b} keys %allsub;
mkdir "$opts{s}", 0755;

open(OUT, ">$opts{s}/$opts{n}_class_raw.tsv") || die "$! $opts{s}/$opts{n}_class_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ac) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $class->{$sub}->{$s} ? $class->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">$opts{s}/$opts{n}_level1_raw.tsv") || die "$! $opts{s}/$opts{n}_level1_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al1) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl1->{$sub}->{$s} ? $lvl1->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">$opts{s}/$opts{n}_level2_raw.tsv") || die "$! $opts{s}/$opts{n}_level2_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al2) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl2->{$sub}->{$s} ? $lvl2->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">$opts{s}/$opts{n}_$opts{s}_raw.tsv") || die "$! $opts{s}/$opts{n}_$opts{s}_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ass) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $ss->{$sub}->{$s} ? $ss->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">$opts{s}/$opts{n}_all_raw.tsv") || die "$! $opts{s}/$opts{n}_all_raw.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$aall) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $all->{$sub}->{$s} ? $all->{$sub}->{$s} : 0;
	}
	print OUT "\n";
}
close OUT;

## NORMALIZED COUNTS

# we divide the totals by 1e6 so that when we do the division we now are normalized per million reads
# NORMALIZED TO ALL READS THAT HAVE A MATCH

map {$total{$_} /= 1e6} keys %total;

open(OUT, ">$opts{s}/$opts{n}_class_norm_all.tsv") || die "$! $opts{s}/$opts{n}_class_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ac) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $class->{$sub}->{$s} ? $class->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">$opts{s}/$opts{n}_level1_norm_all.tsv") || die "$! $opts{s}/$opts{n}_level1_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al1) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl1->{$sub}->{$s} ? $lvl1->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">$opts{s}/$opts{n}_level2_norm_all.tsv") || die "$! $opts{s}/$opts{n}_level2_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al2) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl2->{$sub}->{$s} ? $lvl2->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">$opts{s}/$opts{n}_$opts{s}_norm_all.tsv") || die "$! $opts{s}/$opts{n}_$opts{s}_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ass) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $ss->{$sub}->{$s} ? $ss->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">$opts{s}/$opts{n}_all_norm_all.tsv") || die "$! $opts{s}/$opts{n}_all_norm_all.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$aall) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $all->{$sub}->{$s} ? $all->{$sub}->{$s}/$total{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

# we divide the totals by 1e6 so that when we do the division we now are normalized per million reads
# NORMALIZED TO ALL READS WITH SUBSYSTEMS

map {$sstotal{$_} /= 1e6} keys %sstotal;

open(OUT, ">$opts{s}/$opts{n}_class_norm_ss.tsv") || die "$! $opts{s}/$opts{n}_class_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ac) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $class->{$sub}->{$s} ? $class->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">$opts{s}/$opts{n}_level1_norm_ss.tsv") || die "$! $opts{s}/$opts{n}_level1_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al1) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl1->{$sub}->{$s} ? $lvl1->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;

open(OUT, ">$opts{s}/$opts{n}_level2_norm_ss.tsv") || die "$! $opts{s}/$opts{n}_level2_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$al2) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $lvl2->{$sub}->{$s} ? $lvl2->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">$opts{s}/$opts{n}_$opts{s}_norm_ss.tsv") || die "$! $opts{s}/$opts{n}_$opts{s}_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$ass) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $ss->{$sub}->{$s} ? $ss->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">$opts{s}/$opts{n}_all_norm_ss.tsv") || die "$! $opts{s}/$opts{n}_all_norm_ss.tsv";
print OUT "\t", join("\t", @allsub), "\n";
foreach my $s (sort {$a cmp $b} keys %$aall) {
	print OUT $s;
	foreach my $sub (@allsub) {
		print OUT "\t";
		print OUT $all->{$sub}->{$s} ? $all->{$sub}->{$s}/$sstotal{$sub} : 0;
	}
	print OUT "\n";
}
close OUT;


open(OUT, ">$opts{s}/README.md") || die "$! $opts{s}/README.md";
print OUT <<END;

# NORMALIZATIONS

Currently we perform three normalizations:

1. \*\_raw.tsv

This is the non-normalised data, so just the raw counts. For each sequence, if it appears in one subsystem we incremenet that count by 1, but if it occurs in more than one subsystem, we increment that count by 1/n (1/2 for 2 subsystems, 1/3 for 3 subsystems, etc).

2. \*\_norm_all.tsv

This data is normalised for _all_ reads, regardless of whether they are in a subsystem or not. This makes smaller numbers. 

3. \*\_norm_ss.tsv

This data is normalised only to the number of reads that match to subsystems, so if there is a lot of other stuff we ignore it.


END
close OUT;
