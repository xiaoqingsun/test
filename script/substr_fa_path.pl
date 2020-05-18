#!/usr/bin/perl -w
use strict;

die "use <path>" unless @ARGV==1;

my $path=shift;

my @all_sites;

my $gxyfile="$path/01products/01GXY/GXY_sample_rigth_info.txt";


if(-e $gxyfile){
	open(IN,$gxyfile)||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		next if(/^#/ || /^Flag/);
		next if($t[2] eq '-');
		if($t[4] eq '.'){$t[4]=$t[3];}
		push @all_sites,["gxy",$t[1],$t[2],$t[3],$t[4],$t[8]];
	}
}

my @filterfile=glob("$path/01products/*/*right.txt");

foreach my $file(@filterfile){
	my $pro=(split /\//,$file)[-2];
	open(IN,$file)||die;
	while(<IN>){
                chomp;
                my @t=split /\t/;
		next if(/^#/ || /^Flag/);
                next if($t[2] eq '-');
                if($t[4] eq '.'){$t[4]=$t[3];}
                push @all_sites,[$pro,$t[1],$t[2],$t[3],$t[4],$t[8]];
        }
	close IN;
}

@all_sites=sort{$a->[1] cmp $b->[1]}@all_sites;
open(OUT,">$path/01products/substr_fa_position")||die;
foreach my $p(@all_sites){
	#my $st=$p->[3]-200;
	my $st=$p->[3]-50; ### panqi 181226
	$st=$st>0?$st:0;
	#my $en=$p->[4]+200;
	my $en=$p->[4]+50; ### panqi 181226
	print OUT "$p->[0]\t$p->[1]\t$p->[2]\t$st\t$en\t$p->[5]\t$p->[3]\t$p->[4]\n";
}
close OUT;

`perl  /PUBLIC/pipline/script/Report/XKFW/V3_180824/script/substr_fa.pl   /PUBLIC/pipline/database/ref/hg19_24chr/human_g1k_v37c_decoy.fasta  $path/01products/substr_fa_position  > $path/check/substr_fa_position.fa`;
