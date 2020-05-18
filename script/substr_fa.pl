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
	next if(exists $hash{"$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6] $cut[7]"});
	my $start=$cut[3]-1;
	my $site=$cut[6]-1;
	my $end=$cut[7];
	my $len=$cut[7]-$cut[6]+1;
	my $seq1=substr($seq{$cut[2]},$start,50);
	my $seq2=substr($seq{$cut[2]},$site,$len);
	my $seq3=substr($seq{$cut[2]},$end,50);
	$hash{"$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6] $cut[7]"}=1;
	my $line1 = ">$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6] $cut[7] [+] -->\n$seq1\[$seq2\]$seq3\n";
	my $seq1_2 = reverse($seq1);
	my $seq2_2 = reverse($seq2);
	my $seq3_2 = reverse($seq3);
	my $line3 = ">$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6] $cut[7] [+] <--\n$seq3_2\[$seq2_2\]$seq1_2\n";
	$seq1_2 =~ tr/ACGTNacgtn/TGCANtgcan/;
	$seq2_2 =~ tr/ACGTNacgtn/TGCANtgcan/;
	$seq3_2 =~ tr/ACGTNacgtn/TGCANtgcan/;
	my $line2 = ">$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6] $cut[7] [-] <--\n$seq3_2\[$seq2_2\]$seq1_2\n";
	$seq1 =~ tr/ACGTNacgtn/TGCANtgcan/;
	$seq2 =~ tr/ACGTNacgtn/TGCANtgcan/;
	$seq3 =~ tr/ACGTNacgtn/TGCANtgcan/;
	my $line4 = ">$cut[0] $cut[1] $cut[2] $cut[3] $cut[4] $cut[5] $cut[6] $cut[7] [-] -->\n$seq1\[$seq2\]$seq3\n";
	print $line1.$line2.$line3.$line4;
}
close IN;




