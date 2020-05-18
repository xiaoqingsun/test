use strict;
use warnings;
use utf8;
use Encode;
use Getopt::Long;
use POSIX;
use Data::Dumper;
use lib '/PUBLIC/pipline/script/Report/XKFW/V3_180824/script/';
use Cardio_03_Yaomim_SUMMARY_gai; #zhouwei change 20170713    yangsheng change 20170921
use Cardio_03_Yaomim_REPORT_gai; #zhouwei change 20170713     yangsheng change 20170921

my ($gresult,$yresult,$dresult,$sample,$outdir,$part,$module);
GetOptions(
	'gresult=s'=> \$gresult,
	'yresult=s'=> \$yresult,
	'dresult=s'=> \$dresult,
	'sample=s' => \$sample,
	'outdir=s' => \$outdir,
	'part=s' => \$part,
    'module=s'=> \$module,   #yangsheng change 20170921
);

my $tex1= Cardio_03_Yaomim_SUMMARY::summary($gresult,$yresult,$dresult,$part,$module); #zhouwei change 20170713 yangsheng change 20170921
my $tex2= Cardio_03_Yaomim_REPORT::report($gresult,$yresult,$dresult,$part); #zhouwei change 20170713  
my $tex="$tex1"."$tex2";
open OUT, ">${outdir}/${sample}_Cardio_Yaomim.tex" or die; #yangsheng change 20170921
print OUT "$tex\n";
