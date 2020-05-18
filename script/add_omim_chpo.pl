#!/usr/bin/perl -w
use strict;
use warnings;
use utf8;
use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Getopt::Long;


my ($path,$out);
GetOptions(
        'path=s'  => \$path,
        'out=s' => \$out,
);

my $usage=<<END;

perl $0 
        -path   需要检测样本的路径list文件，
        -out    输出路径,输出相应产品的true的点，且YC，GXY，CS会与最新位点数据库匹配,返回与数据库匹配的结果
END


my  %omim_list;
my  %gene_list;
my  %gene_f;

my $cusigene="/PUBLIC/pipline/database/CHPO/omim.sort.gene2omim2hpo.add_system.txt";
open(IN,$cusigene)||die;
my $head2=<IN>;
while(<IN>){
	chomp;
	my @t=split /\t/;
        my @genes=split /, /,$t[2];
	if($t[1] ne ''){
		push @{$gene_list{$t[1]}},$_;
                $gene_f{$t[1]}=1;
	}
        foreach my $g(@genes){
                next if($g eq $t[1]);
		push @{$gene_list{$g}},$_;
		$gene_f{$g}=1;
	}
	push @{$omim_list{$t[0]}},$_;
        $gene_f{$t[0]}=1;
}
close IN;

open(IN,$path)||die;
my $basename=basename $path;
open(OUT,">$out/$basename.anno")||die;
my $head=<IN>;
chomp $head;
my @a=split /\t/,$head;
my $nhead= @a;
print OUT "$head\t$head2";
while(<IN>){
	chomp;
	my @t=split /\t/;
	my @genes=split /;/,$t[8]; 
        my $flag=0;
        my $nline=@t;
        my $nbu=$nhead-$nline;
        my $bug="\t"x$nbu;
        foreach my $G(@genes){
		last if($flag==1);
		if( exists $gene_f{$G} ){
			$flag=1;
			foreach my $p(@{$gene_list{$G}}){
				print OUT "$_\t$bug$p\n";					
			}
		}
	}
	if($flag==0){
		if( exists $gene_f{$t[66]}){
                        $flag=1;
			foreach my $p(@{$omim_list{$t[66]}}){
                                print OUT "$_\t$bug$p\n";                               
                        }
		}
	}
	if($flag==0){
		print OUT "$_\t$bug","No_anno\n";
	}
}
close IN;
close OUT;
