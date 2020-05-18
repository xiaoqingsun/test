use strict;
use warnings;
use utf8;
use Encode;
use Getopt::Long;
use POSIX;
use Data::Dumper;

my ($sample,$outdir,$part,$yang_gxy,$yang_cs,$yang_dm,$yang_wes,$site_gxy,$site_cs,$site_dm,$site_wes);
GetOptions(
        'sample=s' => \$sample,
        'outdir=s' => \$outdir,
        'part=s' => \$part,
	'yang_gxy=i' => \$yang_gxy,
	'yang_cs=i' => \$yang_cs,
        'yang_dm=i' => \$yang_dm,
        'yang_wes=i' => \$yang_wes,
	'site_gxy=s' => \$site_gxy,
	'site_cs=s' => \$site_cs,
        'site_dm=s' => \$site_dm,
        'site_wes=s' => \$site_wes,
);


open OUT, ">${outdir}/${sample}_Cardio_Sanger.tex" or die; #sunxq change 20171220

my $cover;
my $fig;
my $gxy_site;
my $cs_site;
my $dm_site;
my $wes_site;
if($yang_gxy==1){
	my %hash;
	open(IN,$site_gxy)||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		if($t[0] eq $sample){
			my $aa="$t[85]：$t[86]，$t[87]";
			$aa=~s/_/\\_/g;
			
			if(!exists $hash{$aa})
			{
				$gxy_site.=qq[\\$aa  \\\\ \n];
				$fig="$outdir/fig/$t[0]-$t[1]-$t[2].jpg";
				if(-e $fig)
				{
					$gxy_site.=qq[ \\\\ \\centerline{\\includegraphics[width=15cm,height=4cm]{$outdir/fig/$t[0]-$t[1]-$t[2].jpg}}  \\\\ \\\\ \\\\ \\\\ \n ];
				}
				else{
					$gxy_site.=qq[\\\\{参照序列\\\\\\\\样本序列\\\\} \n \\\\ \n \\downarrow \n ];
				}
			}
			$hash{$aa}=1;
		}
	}
	close IN;
}
if($yang_cs==1){
	open(IN,$site_cs)||die;
	my %hash;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		if($t[0] eq $sample){
			my $aa="$t[85]：$t[86]，$t[87]";
			$aa=~s/_/\\_/g;
			
			if(!exists $hash{$aa})
			{
				$cs_site.=qq[\\$aa  \\\\ \n];
				$fig="$outdir/fig/$t[0]-$t[1]-$t[2].jpg";
				if(-e $fig)
				{
					$cs_site.=qq[ \\\\ \\centerline{\\includegraphics[width=15cm,height=4cm]{$outdir/fig/$t[0]-$t[1]-$t[2].jpg}}  \\\\ \\\\ \\\\ \\\\ \n ];
				}
				else{
					$cs_site.=qq[\\\\{参照序列\\\\\\\\样本序列\\\\} \n \\\\ \n \\downarrow \n ];
				}
			}
			$hash{$aa}=1;
		}
	}
	close IN;
}

if($yang_dm==1){
	open(IN,$site_dm)||die;
	my %hash;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		if($t[0] eq $sample){
			my $aa="$t[85]：$t[86]，$t[87]";
			$aa=~s/_/\\_/g;
			
			if(!exists $hash{$aa})
			{
				$dm_site.=qq[\\$aa  \\\\ \n];
				$fig="$outdir/fig/$t[0]-$t[1]-$t[2].jpg";
				if(-e $fig)
				{
					$dm_site.=qq[ \\\\ \\centerline{\\includegraphics[width=15cm,height=4cm]{$outdir/fig/$t[0]-$t[1]-$t[2].jpg}}  \\\\ \\\\ \\\\ \\\\ \n ];
				}
				else{
					$dm_site.=qq[\\\\{参照序列\\\\\\\\样本序列\\\\} \n \\\\ \n \\downarrow \n ];
				}
			}
			$hash{$aa}=1;
		}
	}
	close IN;
}

if($yang_wes==1){
	open(IN,$site_wes)||die;
	my %hash;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		if($t[0] eq $sample){
			my $aa="$t[85]：$t[86]，$t[87]";
			$aa=~s/_/\\_/g;
			
			if(!exists $hash{$aa})
			{
				$wes_site.=qq[\\$aa  \\\\ \n];
				$fig="$outdir/fig/$t[0]-$t[1]-$t[2].jpg";
				if(-e $fig)
				{
					$wes_site.=qq[ \\\\ \\centerline{\\includegraphics[width=15cm,height=4cm]{$outdir/fig/$t[0]-$t[1]-$t[2].jpg}}  \\\\ \\\\ \\\\ \\\\ \n ];
				}
				else{
					$wes_site.=qq[\\\\{参照序列\\\\\\\\样本序列\\\\} \n \\\\ \n \\downarrow \n ];
				}
			}
			$hash{$aa}=1;
		}
	}
	close IN;
}

if($yang_gxy==1 || $yang_cs==1 || $yang_dm==1  ||  $yang_wes==1)
{ 
        my $flag=0;
	$cover=qq[%%%%%%%%%%%%%%%%% sanger验证 结果 %%%%%%%%%%%%%%%
\\newpage
\\addcontentsline{toc}{section}{Part$part — 位点一代验证结果}
\\ThisTileWallPaper{\\paperwidth}{\\paperheight}{XKhead.pdf}
\\begin{center}
\\zihao{2}{\\color{DarkBlue} \\sye \\sym Part$part — 位点一代验证结果} \\
\\end{center}
\\par \\vspace{1cm}
];
       if($yang_gxy==1){
		$cover.=qq[\\includegraphics[height=0.4cm]{XKmini.pdf} \\zihao{4}{\\color{DarkBlue} \\sym{ 单基因高血压/低血钾位点}}\\par
\\fontsize{10}{15}\\selectfont
\\$gxy_site
%\\\\\\\\{参照序列}\\\\\\\\{样本序列}\\\\
%\\fontsize{24}{25}
%\\\\
%\\\\\\downarrow
];
                $flag++;
       }
       if($yang_cs==1){
                if($flag>0){$cover.=qq[\\newpage
];}
                $cover.=qq[\\includegraphics[height=0.4cm]{XKmini.pdf} \\zihao{4}{\\color{DarkBlue} \\sym{ 遗传性心血管病位点}}\\par
\\fontsize{10}{15}\\selectfont
\\$cs_site
%\\\\\\\\{参照序列}\\\\\\\\{样本序列}\\\\
%\\fontsize{24}{25}
%\\\\
%\\\\\\downarrow
\\par
];
		$flag++;
       }
       if($yang_dm==1){
                if($flag>0){$cover.=qq[\\newpage
];}
                $cover.=qq[\\includegraphics[height=0.4cm]{XKmini.pdf} \\zihao{4}{\\color{DarkBlue} \\sym{ 单基因糖尿病位点}}\\par
\\fontsize{10}{15}\\selectfont
\\$dm_site
%\\\\\\\\{参照序列}\\\\\\\\{样本序列}\\\\
%\\fontsize{24}{25}
%\\\\
%\\\\\\downarrow
\\par
];
                $flag++;
       }
	if($yang_wes==1){
                if($flag>0){$cover.=qq[\\newpage
];}
                $cover.=qq[\\includegraphics[height=0.4cm]{XKmini.pdf} \\zihao{4}{\\color{DarkBlue} \\sym{ 全外显子检测位点}}\\par
\\fontsize{10}{15}\\selectfont
\\$wes_site
%\\\\\\\\{参照序列}\\\\\\\\{样本序列}\\\\
%\\fontsize{24}{25}
%\\\\
%\\\\\\downarrow
\\par
];
                $flag++;
       }
}

print OUT "$cover\n";

close OUT;
