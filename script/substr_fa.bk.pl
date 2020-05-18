#! usr/bin/perl -w
use strict;

die "use<fa><regin_list>\n" unless @ARGV==2;
my %seq;
open(IN,$ARGV[0])||die;
$/ = ">";
<IN>;
while(<IN>)
{
	/(.+)\n/;
	my $id = (split /\s+/)[0];
	s/.+\n//;
	s/\s+|>//g;
	$seq{$id}=$_;
	
}
close IN;

my %hash;

$/ = "\n";
open(IN,$ARGV[1])||die;
while(<IN>)
{
	chomp;
	my @cut=split /\s+/;
	next if(exists $hash{"$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6]"});
	my $start=$cut[3]-1;
	my $end=$cut[4]-1;
	my $len=$cut[4]-$cut[3]+1;
	my $seq=substr($seq{$cut[2]},$start,$len);
	$hash{"$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6]"}=1;
	print ">$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6]\n$seq\n";
}
close IN;




