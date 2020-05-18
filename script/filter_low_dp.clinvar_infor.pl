my $file1 = shift;
my $file2 = shift;
my $file3 = shift;

my @cut;
my %tmp;
my $line="";
my $len=0;
open IN, $file2 or die $!;
while(<IN>)
{
	chomp;
	@cut=split /\t/,$_;
	$tmp{$cut[0]}=$_;
}
close IN;

open IN, $file3 or die $!;
while(<IN>)
{
	chomp;
	@cut=split /\t/,$_;
	$tmp{$cut[0]}=$_;
}
close IN;

$line=$tmp{"#AlleleID"};
print "CHR\tStart\tEnd\tLength\tdepth\tchr2\tstart2\tend2\tGENE\tmRNA\t+/-\tEXON\tCDS\tEXON_ALL\texon_length\t$line\n";
open IN, $file1 or die $!;
while(<IN>){
	chomp;
	@cut=split /\t/,$_;

	if($cut[19] eq "-1")
	{
		print $_."\n";
	}
	else{
		next if($cut[19]!~/^Likely pathogenic/ && $cut[19]!~/^Pathogenic/);
		$line=join("\t",@cut[0..14]);
		print $line."\t$tmp{$cut[18]}\n";
	}
}
close IN;
