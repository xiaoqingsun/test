#!/usr/bin/perl -w
use strict;
use warnings;
use utf8;
use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Encode;

my ($out,$types,$reference);
GetOptions(
        'out=s' => \$out,
        'type=s' => \$types,
        'reference=s' => \$reference,
);

my (%allref,$reportfile,$finalfile,$refcol);

open(IN,$reference)||die;
while(<IN>){
	chomp;
	my @t=split /\t/;
	$allref{$t[0]}=$t[1];
}
close IN;

$refcol=82;
if($types eq 'GXY'){
	$reportfile="$out/GXY_sample_rigth_to_PDF_final.txt";
}else{
	$reportfile="$out/$types.pdf.txt";
}

$finalfile="$reportfile.ref";

open(IN,"<:utf8",$reportfile)||die;
open(OUT,">$finalfile")||die;
my $count=0;
while(<IN>){
	chomp;
	my @t=split /\t/;
	$count++;
	my %refs;
	if($count==1){
		if($t[$refcol] ne '参考文献'){
			die"请确认参考文献的列数\t$refcol\t$t[$refcol]\n";		
		}
		print OUT "$_\n";
	}else{
		if($t[$refcol] eq '-' ||$t[$refcol] eq '.'){
			print OUT "$_\n";
			next;
		}
		my @a=split /\|\|/,$t[$refcol];
		#my @a1=split /\|/,$a[0];
		#my @a2=split /\|/,$a[1];
		foreach my $aa(@a){
			my @a1=split /\|/,$aa;
			foreach my $r(@a1){
				if($r=~/^\d+$/){
					if(!exists $allref{$r}){
						print "请确认参考文献pubmedID:$r\n";
					}else{
						my $subref=$allref{$r};
						$refs{pub}{$subref}=1;
					}
				}else{
					$refs{ref}{$r}=1;
				}
			}
		}
		my @references=();
		foreach my $r1(keys %{$refs{pub}}){
			push @references,$r1;
		}
		foreach my $r1(keys %{$refs{ref}}){
                        push @references,$r1;
                }
		my $n=@references;
		my $ref_line='-';
		$ref_line=join("|",@references) if($n>0);
		$t[$refcol]=$ref_line;
		my $line=join("\t",@t);
		print OUT "$line\n";
	}
}
close IN;
close OUT;
`mv $reportfile $reportfile.bak`;
`mv $finalfile $reportfile`;
