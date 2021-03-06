use strict;
use warnings;
#use utf8;
#use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Smart::Comments;
use Getopt::Long;
use Spreadsheet::XLSX;
use List::Util qw/max min sum maxstr minstr shuffle/;


my ($samples,$site_file,$out,$types,$qc_stat,$hospital,$usage,$genelist)=('','','','','','','','');
GetOptions(
    'sample=s'=> \$samples,
    'hospital=s'=> \$hospital,
    'site=s' => \$site_file,
    'out=s'  => \$out,
    'type=s' => \$types,
    'qc=s'   => \$qc_stat,
    'genelist=s'=> \$genelist,  #sxq 5.29
	#'help|h' => \$usage,
);


$usage=<<END;

perl $0
	-sample 样本ID
	-site   需要检测样本的位点注释文件
	-type   产品类型, GXY/CS/DM/WES
	-qc     样本质控文件
	-out    样本报告输出路径
	-genelist 检测疾病基因列表
END

if($samples eq '' or $site_file eq '' or $types eq '' or $qc_stat eq '')
{
	die "$usage\n";
}

chomp($out=`pwd`) if($out eq '');
$out =~ s/\/$//;

my @cut;
my ($name,$sex,$year,$id_len,$shizi,$clinical_info,$jiazu,$yiyuan,$header,$footer)=('','','',0,"口腔拭子","无","无","NULL",'','');
my ($yemei,$biaoti,$biaoti2,$info_sample,$info1,$part1,$part2,$part3,$part4,$part5,$part6,$part6_0,$part7_1,$part7,$part8,$check);
my $part1_list='';
my $negative=0;
my $negative2=0;
my $fb_num=0;
#my $bin = "/PUBLIC/work/panqi/diagnosis_report/script";
my $bin = "/PUBLIC/pipline/script/Report/XKFW/V2_171130/script";


################### part 0 sample info #####################
my $hutongbiao="/PUBLIC/pipline/database/sheet_hutong/";
my ($file,$file2) = judge_lastest_sitedb($hutongbiao);
`perl $bin/spectrum_XLSX.pl $file` unless(-e $file2);

#open IN, "<:encoding(utf-8)", "$file2" or die $!;
open IN, "$file2" or die $!;
<IN>;
while(<IN>)
{
	chomp;
	#next if(/^序号/)；
	@cut=split /\t/, $_;
	if($samples eq $cut[2])
	{
		$name=$cut[4];
		$sex=$cut[7];
		$year=$cut[8];
        if($year=~/岁/ or $year=~/月/){print $year;}
        else{$year=$year."岁"}#wph 180823
		$yiyuan=$cut[9];
		$cut[1]=~s/\//-/g;
		
		$id_len=length $cut[3];
		if($id_len > 8){$shizi="外周血";}
		if($shizi eq "口腔拭子")
		{
			$info1=qq[\\ \\ \\rlap{\${\\surd}\$}\${\\square}\$ 口腔拭子 \\ \\ \${\\square}\$ 外周血];
		}
		elsif($shizi eq "外周血")
		{
			$info1=qq[\\ \\ \${\\square}\$ 口腔拭子 \\ \\ \\rlap{\${\\surd}\$}\${\\square}\$ 外周血];
		}

		$info_sample=
qq[{\\sym{姓名：}} $name & {\\sym{性别：}} $sex & {\\sym{年龄：}} $year & {\\sym{收样日期：}} $cut[1] \\\\ \\hline
\\multicolumn{2}{l} {{\\sym{样本类型：}} $info1} & \\multicolumn{2}{l} {{\\sym{检测方法：}}高通量测序} \\\\ \\hline
\\end{tabular}
];	
		#\\multicolumn{2}{l} {{\\sym{临床诊断：}} $clinical_info } & \\multicolumn{2}{l} {{\\sym{家族史：}} $jiazu } \\\\ \\hline	
	}
}
close IN; 

open FILE, $hospital or die $!;
while(<FILE>)
{
	next if(/^product/);
	chomp;
	@cut=split /\t/,$_;
	if($yiyuan eq "糖友" &&  $yiyuan eq $cut[0])
	{
		$header=$cut[1]."\n".$cut[2]."\n";
		$footer=$cut[3]."\n";
	}
	elsif($cut[0] eq "诺禾心康")
	{
		$header=$cut[1]."\n".$cut[2]."\n";
		$footer=$cut[3]."\n";
	}
}
close FILE;

################### part 1 site info #################

my %het_hom = ('het'=>'杂合','hom'=>'纯合','unknown'=>'未明');
my %level_value = ('致病'=>1,'可能致病'=>2,'临床意义未明1级'=>3,'临床意义未明2级'=>4,'临床意义未明3级'=>5,'临床意义未明4级'=>6,'临床意义未明5级'=>7,'可能良性'=>8,'良性'=>9);
my %mute_des = 
   ('frameshift deletion'=>'移码缺失变异',
    'frameshift insertion'=>'移码插入变异',
    'frameshift substitution'=>'移码替换变异',
    'nonframeshift deletion'=>'非移码缺失变异',
    'nonframeshift insertion'=>'非移码插入变异',
    'nonframeshift substitution'=>'非移码替换变异',
    'stopgain'=>'无义变异',
    'stopgain2'=>'缺失移码变异',
    'stoploss SNV'=>'终止缺失变异',
    'synonymous SNV'=>'同义变异',
    'nonsynonymous SNV'=>'错义变异',
    'SNV'=>'单核苷酸变异',
    'splicing SNV'=>'剪切区变异');

#   'splicing SNV'=>'剪切区单核苷酸变异');

my %index;
my (%gene_site,%gene_site_level,%gene_site_level_fb,%disease_site_level,%gene_level,%gene_level_fb,%gene_site_level2,%gene_level2,%gene_level2_fb);
my ($disease_name,$disease_gene);
my (%site_dis,%site_dis_bz,%site_dis_fb);
my ($site_info,$exon);
my ($c1,$c2,$c3,$c4,$c5,$c6,$c6_2,$c7,$c8);
my ($sn,$s1,$s2,$s3)=(0,'','','');
my (@paper,@arr,@arr2);
my $num_line=0;
my $num_line_yang=0;
my $person="NULL";
my %second_line=();
my %site_info_tmp=();

#open IN, "<:utf8", "$site_file" or die $!;

open IN, $site_file or die $!;
while(<IN>)
{
	chomp;
	$_=~s/\_/\\\_/g;
	$_=~s/\%/\\\%/g;
	@cut=split /\t/,$_;

	if($num_line == 0)
	{
		for(my $i=0;$i<@cut;$i++)
		{
			$index{$cut[$i]}=$i;
		}		
	}
	else{
		next unless ($cut[0]=~/$samples/);
		next if($cut[0]=~/^####/);
		next if($cut[0]=~/^NN/);
		$person=$cut[$index{'解读人'}] if(exists $index{'解读人'} && $cut[$index{'解读人'}] ne "-");
		$person="NULL" if(!exists $index{'解读人'});
		$site_info=join("\_",@cut[1..3,5,6]);

		if($site_info eq "-_-_-_-_-")
		{
			$negative=2;
			next;
		}

		if($cut[$index{'PDF_NMID'}] eq "\.\/\/")
		{
			$exon="NULL";
		}else{
			$exon=join "\\\\", ( split /:/, $cut[$index{'PDF_NMID'}] )[0,1];
		}
		if(!exists $gene_site{$cut[7]}{$site_info})
		{
			$gene_site{$cut[7]}{$site_info}="$cut[$index{'HGNC基因名'}]\t$cut[$index{'CDS'}]\t$cut[$index{'PEP'}]\t$exon\t$het_hom{$cut[31]}\t$mute_des{$cut[9]}\t$cut[21]\t$cut[$index{'基因描述'}]\t$cut[$index{'位点描述'}]" if(exists $index{'HGNC基因名'});
			$gene_site{$cut[7]}{$site_info}="$cut[$index{'基因名称'}]\t$cut[$index{'CDS'}]\t$cut[$index{'PEP'}]\t$exon\t$het_hom{$cut[31]}\t$mute_des{$cut[9]}\t$cut[21]\t$cut[$index{'基因描述'}]\t$cut[$index{'位点描述'}]" if(exists $index{'基因名称'});
		}

		if($cut[0]=~/^###/)
		{
			push @{ $site_dis_bz{$cut[7]}{$site_info} }, [ $cut[$index{'疾病名称'}], $cut[$index{'遗传模式'}], $cut[$index{'ACMG条目'}], $cut[$index{'致病性结论'}] ]  if(exists $index{'ACMG条目'});
			push @{ $site_dis_bz{$cut[7]}{$site_info} }, [ $cut[$index{'疾病名称'}], $cut[$index{'遗传模式'}], $cut[$index{'ACMG对应条目'}], $cut[$index{'位点结论'}] ]  if(exists $index{'ACMG对应条目'});
		}
		elsif($cut[0]=~/^##/)
		{
			$negative2=1;
			if(exists $index{'ACMG条目'})
			{
				$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{ $cut[$index{'致病性结论'}] };
				push @{$gene_level_fb{$cut[7]}},$level_value{ $cut[$index{'致病性结论'}] };
				push @{ $site_dis_fb{$cut[7]}{$site_info} }, [ $cut[$index{'疾病名称'}], $cut[$index{'遗传模式'}], $cut[$index{'ACMG条目'}], $cut[$index{'致病性结论'}] ];
			}
			if(exists $index{'ACMG对应条目'})
			{
				$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{$cut[$index{'位点结论'}]};
				push @{$gene_level_fb{$cut[7]}},$level_value{$cut[$index{'位点结论'}]};
				push @{ $site_dis_fb{$cut[7]}{$site_info} }, [ $cut[$index{'疾病名称'}], $cut[$index{'遗传模式'}], $cut[$index{'ACMG对应条目'}], $cut[$index{'位点结论'}] ];
			}
		}
		else{
			if(exists $index{'位点结论'})
			{
				push @{ $gene_site_level{$cut[7]}{$site_info} }, $level_value{ $cut[$index{'位点结论'}]  };
				push @{ $site_dis{$cut[7]}{$site_info} }, [ $cut[$index{'位点结论'}], $cut[$index{'疾病名称'}], $cut[$index{'遗传模式'}], $cut[$index{'疾病描述'}] ];
				$disease_site_level{$cut[7]}{$site_info}{$cut[$index{'疾病名称'}]}=$cut[$index{'位点结论'}];
			}
			elsif(exists $index{'致病性结论'})
			{
				push @{ $gene_site_level{$cut[7]}{$site_info} }, $level_value{ $cut[$index{'致病性结论'}]  };
				push @{ $site_dis{$cut[7]}{$site_info} }, [ $cut[$index{'致病性结论'}], $cut[$index{'疾病名称'}], $cut[$index{'遗传模式'}], $cut[$index{'疾病描述'}] ];
				$disease_site_level{$cut[7]}{$site_info}{$cut[$index{'疾病名称'}]}=$cut[$index{'致病性结论'}];
			}
			if($cut[$index{'参考文献'}] ne '-' && $cut[$index{'参考文献'}] ne '.')
			{ 
				push @paper, split(/\|/, $cut[$index{'参考文献'}]); 
			}
			if(!exists $site_info_tmp{$site_info})
			{
				$num_line_yang++;
				$site_info_tmp{$site_info}=1;
			}
			#$person=$cut[$index{'解读人'}] if(exists $index{'解读人'});
			#$person="NULL" if(!exists $index{'解读人'});
		}
	}
	$num_line++;
}
close IN;

################### paper #################
my %arr_tmp;
my $part8_1='';
$part8='';
if(@paper)
{
	for(my $i=0;$i<@paper;$i++)
	{
		if(!exists $arr_tmp{$paper[$i]})
		{
			$part8_1.="\\item $paper[$i]" if($paper[$i] ne '-' && $paper[$i] ne '.');
			$arr_tmp{$paper[$i]} = 1;
		}
	}

	$part8=<<'END';
\color{MyFontGray}
\vspace*{-3mm}
\zihao{6}{{\sym 参考文献}} \\
\begin{spacing}{1.4}
\zihao{6} \color{MyFontGray}
\setdefaultleftmargin{2em}{}{}{}{.5em}{.5em}
\begin{compactenum}
END

	$part8.=$part8_1."\n\\end{compactenum} \n \\end{spacing} \n";
}

################### site #################
$part3='';
if(%site_dis_bz)
{
$part3='{\\sym{备注：}} ';
foreach my $key1 (sort {$a cmp $b} keys %site_dis_bz)
{
	$part3.="$key1基因可能还与";
	foreach my $key2 (sort {$a cmp $b} keys %{$site_dis_bz{$key1}})
	{
		@arr=@{$site_dis_bz{$key1}{$key2}};
		push @arr2, @arr;
	}
	
	for(my $i=0;$i<@arr2;$i++)
	{
		if(!exists $arr_tmp{$key1}{$arr2[$i][0]})
		{
			$part3.="$arr2[$i][0]($arr2[$i][1])" if($i == 0);
			$part3.="、$arr2[$i][0]($arr2[$i][1])" if($i > 0);
			$arr_tmp{$key1}{$arr2[$i][0]}=1;
		}
	}
	$part3.="疾病相关，但目前研究较少，如已经出现该病的临床表型，请结合临床考虑该位点的致病性。";
}
$part3.="\\\\ \\\\ \n";
}

if(%gene_site_level)
{

$part1=<<'END';
\noindent \\ \\ \\
{\sym{检测结果：}} \\
\renewcommand\arraystretch{1.5}
END
$part1.=qq[\\fontsize{10.5}{12.6}\\selectfont{\\color{MyFontGray} 受检者检测到 $num_line_yang 个变异位点，具体信息如下表所示：}\\\\];
$part1.=<<'END';
\fontsize{9}{10.8}\selectfont
\begin{longtable}{| C{15mm} | C{16mm} | C{18mm} | C{17mm} | C{8mm} | C{8mm} | C{18mm} | L{35mm} | }
\hline
{\sym{变异基因}} & {\sym{核苷酸变异}} & {\sym{氨基酸变异}} & \makecell[c]{ \sym{转录本}\\ \sym{外显子} } & {\sym{变异状态}} & {\sym{变异类型}} & {\sym{致病性}} & \makecell[c]{\sym{相关疾病(遗传模式)}} \\ \hline
END

foreach my $key1 (sort {$a cmp $b} keys %gene_site_level)
{
	foreach my $key2 (sort {$a cmp $b} keys %{$gene_site_level{$key1}})
	{
		my @level_1=@{ $gene_site_level{$key1}{$key2} };
		my $level_min = min(@level_1);
	
		$gene_site_level2{$key1}{$key2}=$level_min;
		push @{$gene_level{$key1}}, $level_min;
	}
}

foreach my $key1 (sort {$a cmp $b} keys %gene_level)
{
	my @level_1=@{ $gene_level{$key1} };
	my $level_min = min(@level_1);
	$gene_level2{$key1}=$level_min;
}

foreach my $key1 (sort { $gene_level2{$a} <=> $gene_level2{$b} } keys %gene_level2)
{
	foreach my $key2 (sort { $gene_site_level2{$key1}->{$a} <=> $gene_site_level2{$key1}->{$b} } keys %{$gene_site_level2{$key1}})
	{
		($c1,$c2,$c3,$c4,$c5,$c6,$c6_2,$s2,$s3)=split /\t/, $gene_site{$key1}{$key2};
		$sn++;
		foreach my $key3 (sort { $disease_site_level{$key1}{$key2}->{$a} <=> $disease_site_level{$key1}{$key2}->{$b} } keys %{ $disease_site_level{$key1}{$key2} })
		{
			my $j=0;
			my @arr=@{ $site_dis{$key1}{$key2} };
			$c8='';
			for(my $i=0;$i<@arr;$i++)
			{
				next if($arr[$i][1] ne $key3);
				$c7=$arr[$i][0];
				if($arr[$i][2] eq "-")
				{
					$arr[$i][2]="不明";
				}
				$c8.="$arr[$i][1]($arr[$i][2])";
				if($j==0)
				{
					$s1=
qq[{\\sym {位点$sn：}} \\\\
{\\sym{疾病描述：}} \\\\
$arr[$i][1]：$arr[$i][3] \\\\ \\\\ ];
				}
				elsif($j>0)
				{
					$s1.=qq[$arr[$i][1]：$arr[$i][3] \\\\ \\\\ ];
				}
				$part2.=qq[ $c1 & $c2 & $c3 & \\makecell[c]{$c4} & $c5 & $c6 & $c7 & $c8 \\\\ \\hline\n ];
			}
			
		}
		
		$part4.=$s1."\n{\\sym{基因描述：}} ".$s2."\\\\ \\\\ \n{\\sym{位点描述：}} ".$s3."\\\\ \\\\ \n";
		$part5.=
qq[{\\sym{Sanger验证：}} \\\\
$c1:  $c2, $c3 \\\\
参照序列 \\\\
样本序列 \\\\
\\vspace*{4.5cm}
\\\\];
	}
}

}
else{
	$negative=1 if($negative2==1);
}


################### fb #################
my ($fb1,$fb2);
$part6_0='';
if(%site_dis_fb)
{
	$fb_num++ if($negative != 2);
	foreach my $key1 (sort {$a cmp $b} keys %gene_level_fb)
	{
		my @level_1=@{ $gene_level_fb{$key1} };
		my $level_min = min(@level_1);
		$gene_level2_fb{$key1}=$level_min;
	}
	
	foreach my $key1 (sort { $gene_level2_fb{$a} <=> $gene_level2_fb{$b} } keys %gene_level2_fb)
	{
		foreach my $key2 (sort {$gene_site_level_fb{$key1}->{$a} <=> $gene_site_level_fb{$key1}->{$b}} keys %{$gene_site_level_fb{$key1}})
		{
			@arr=split(/\t/, $gene_site{$key1}{$key2});
			@arr2=@{$site_dis_fb{$key1}{$key2}};
			my $i=0;
			my ($acmg,$zbx);
			for($i=0;$i<@arr2;$i++)
			{
				$arr2[$i][0]=~s/\(/\\\\\(/g;
				$arr2[$i][1]="不明" if($arr2[$i][1] eq "-");
				$fb2="$arr2[$i][0]($arr2[$i][1])" if($i==0);
				$fb2.="\\\\$arr2[$i][0]($arr2[$i][1])" if($i>0);#wph del 181010
				$acmg=$arr2[$i][2];
				$acmg=~s/,/, /g;
				$zbx =$arr2[$i][3];
			}
			$fb1.=qq[ $arr[0] & $arr[1] & $arr[2] & \\makecell*[c]{$arr[3]} & $arr[4] & $arr[5] & $acmg & $zbx & \\makecell*[l]{$fb2} \\\\ \\hline\n ];
		}
	}
	$part6_0="{\\sym{附表$fb_num：}} \\\\ \n";
	$part6_0.=<<'END';
下表为无ACMG证据、证据不足、证据矛盾（偏致病类证据与偏良性证据共存）或偏良性证据的位点（可能与相关疾病出现的表型有关，具体情况请结合临床综合分析）\\
\fontsize{7}{8.4}\selectfont
\begin{longtable}{| C{14mm} | C{15mm} | C{18mm} | C{15mm} | C{6mm} | C{6mm} | C{10mm} | C{12mm} | L{40mm} | } \hline
{\sym{变异基因}} & {\sym{核苷酸变异}} & {\sym{氨基酸变异}} & \makecell[c]{ \sym{转录本}\\ \sym{外显子} } & {\sym{变异状态}} & {\sym{变异类型}} & {\sym{ACMG条目}} & {\sym{致病性}} & \makecell[c]{\sym{相关疾病(遗传模式)}} \\ \hline
END

	$part6_0.=$fb1."\\end{longtable} \n";
    print $part6_0;
}

################### IDT_stat #################
if($types eq "WES")
{
	open IN, "$genelist" or die $!;
	my $i=0;
	my @cut_gene;
	my %tmp_gene;
	while(<IN>)
	{
		chomp;
		@cut=split /\t/,$_;
		if($cut[0] eq $samples )
		{
			$disease_name=$cut[1] if($i==0);
			$disease_name.=",$cut[1]" if($i>0);
			@cut_gene=split /, /,$cut[2];
			foreach my $key (@cut_gene)
			{
				$tmp_gene{$key}=1;
			}
			$i++;
		}
	}
	$disease_gene=keys %tmp_gene;
	close IN;
}

################### negative #################
if($negative==2)
{
	$part1=<<'END';
\noindent \\ \\ \\
\begin{spacing}{1.3}
{\sym{检测结果：}} \\
END

	if($types eq "DM")
	{
		$part1.="此次单基因糖尿病基因检测，包含38个疾病相关基因，筛查3类疾病。";
	}
	if($types eq "CS")
	{
		$part1.="此次遗传性心血管病基因检测，包含179个疾病相关基因，筛查36种疾病亚型。";
	}
	if($types eq "GXY")
	{
		$part1.="此次单基因高血压/血钾异常基因检测，包含41个疾病相关基因，筛查15种疾病亚型。";
	}
	if($types eq "WES")
	{
		$part1.="此次全外显子基因检测，包含".$disease_name."疾病相关的".$disease_gene."个基因。";
	}
	$part1.=<<'END';
检测结果显示您并未发生有风险的位点变异，因此，您的疾病症状并不是由这些基因的变异引起的，可能是目前研究尚未涉及的基因变异，也可能是由于其他因素，如不良生活习惯、环境、其他疾病等引起。建议您注意生活细节，积极改变不良生活方式，并配合医生治疗。\\ \\
以上检测结果仅基于当前的科学研究，随研究文献的更新，疾病-基因的对应关系可能会发生变化，因此我们的报告具有一定的时效性。\\
\end{spacing}
\noindent \\ \\
END

}
	
if($negative==1)
{
	$part1=<<'END';
\noindent \\ \\ \\
\begin{spacing}{1.3}
{\sym{检测结果：}} \\
END
	
	if($types eq "DM")
	{
		$part1.="此次单基因糖尿病基因检测，包含38个疾病相关基因，筛查3类疾病。";
	}
	if($types eq "CS")
	{
		$part1.="此次遗传性心血管病基因检测，包含179个疾病相关基因，筛查36种疾病亚型。";
	}
	if($types eq "GXY")
	{
		$part1.="此次单基因高血压/血钾异常基因检测，包含41个疾病相关基因，筛查15种疾病亚型。";
	}
	if($types eq "WES")
	{
		$part1.="此次全外显子基因检测，包含".$disease_name."疾病相关的".$disease_gene."个基因。";
	}
	$part1.=<<'END';
检测结果显示您并未发生高风险的位点变异；一些变异位点由于相关研究较少、ACMG证据条目不足，基于现有研究判断其致病性等级较低，但不排除导致疾病的可能性，这类位点的致病性等级及相关证据条目展示在附表中，供您和临床医生参考。\\ \\
若检出位点的相关疾病与您现有的临床表型不符，那么您的相关临床症状可能与其他未知基因变异、疾病的不完全外显性、不良生活习惯或环境等因素相关。建议您注意生活细节，积极改变不良生活方式，或去医院进行相关疾病的进一步检测，听从医生的临床建议，定期体检，并采取积极的生活干预。\\ \\
若您未发现任何相关疾病表型，则这类 ACMG 致病类证据较少的位点可能并不致病。\\ \\
以上检测结果仅基于当前的科学研究，随研究文献的更新，位点致病性情况可能会发生变化，属于正常现象。如果您想进一步了解基因变异的信息或家族遗传史，建议联系家人一起进行检测。\\ \\
\end{spacing}
\noindent \\ \\
END
	
}

################### list #################
if($negative==1 or $negative==0 or $negative==2)
{
	$fb_num++;
	$part1_list="\\fontsize{10.5}{12.6}\\selectfont\n";
	if($types eq "DM")
	{
		$part1_list.="{\\sym{附表$fb_num：}}";
		$part1_list.=<<'END';
本次检测包含的单基因糖尿病相关基因如下：\\
\renewcommand\arraystretch{1.8}
\fontsize{10.5}{12.6}\selectfont
\begin{longtable}{ L{8.2cm}  L{7.3cm} }
\hline
\makecell[c]{\sym{检测疾病}} & \makecell[c]{\sym{相关基因}} \\ \hline
\makecell[{{p{8.1cm}}}]{年青的成年发病型糖尿病\\MODY(Maturity Onset Diabetes of Young)} & {HNF4A 、GCK 、HNF1A 、PDX1 、HNF1B 、NEUROD1 、KLF11 、CEL 、PAX4 、INS 、BLK 、ABCC8 、KCNJ11 、APPL1 、PCBD1}  \\ \hline
\makecell[{{p{8.1cm}}}]{新生儿糖尿病(NDM)\\Neonatal Diabetes Mellitus} & {KCNJ11 、INS 、GCK 、PDX1 、MNX1 、GLIS3 、IER31P1 、NEUROD1 、NEUROG3 、PAX6 、PTF1A 、RFX6 、GATA4 、NKX2-2 、SLC19A2 、IER3IP1 、PLAGL1 、ZFP57 、ABCC8 、HNF1B 、PLAGL1} \\ \hline
\makecell[{{p{8.1cm}}}]{其他单基因糖尿病} & {SLC2A2 、WFS1 、TRMT10A/PPP1R15B 、EIF2AK3 、CISD2 、INSR 、ALMS1 、BLM 、HADH 、FOXP3} \\ \hline
\end{longtable}
\noindent \\
END
	}
	
	if($types eq "GXY")
	{
		$part1_list.="{\\sym{附表$fb_num：}}";
		$part1_list.=<<'END';
本检测共筛查15种单基因高血压/血钾异常疾病，包含41个相关基因，具体疾病名称请参见下表。 \\
\renewcommand\arraystretch{1.8}
\fontsize{10.5}{12.6}\selectfont
\begin{longtable}{ L{8.3cm}  L{7.2cm} }
\hline
\makecell[c]{\sym{检测疾病}} & \makecell[c]{\sym{相关基因}} \\ \hline
END

		open IN, "/PUBLIC/pipline/script/Report/XKFW/V2_171130/database/gxy_disease_name2.txt" or die $!;
		while(<IN>)
		{
			chomp;
			@cut=split /\t/,$_;
			$part1_list.="\\makecell[{{p{8.2cm}}}]{$cut[0]} & \\makecell[{{p{7cm}}}]{$cut[1]} \\\\ \\hline \n";

			#for(my $i=0;$i<@cut;$i++)
			#{
			#	$part1_list.="\\makecell[{{p{8.2cm}}}]{$cut[0]} & \\makecell[{{p{7cm}}}]{$cut[1]} \\\\ \\hline \n";
			#}
		}
		close IN;
		
		$part1_list.=<<'END';
\end{longtable}
\noindent \\
END

	}
	
	if($types eq "CS")
	{
		$part1_list.="{\\sym{附表$fb_num：}}";
		$part1_list.=<<'END';
本次遗传性心血管病基因检测包括 179 个基因的全部外显子区域，共涉及5大类、 36 种遗传性心血管病，具体疾病名称请参见下表。\\
\renewcommand\arraystretch{1.8}
\fontsize{9}{10.8}\selectfont
\begin{longtable}{ L{8.3cm}  L{7.2cm} }
\hline
END
		
	$part1_list.=read_CS('xinlv');
    	$part1_list.=read_CS('xinji');
    	$part1_list.=read_CS('zhudongmai');
    	$part1_list.=read_CS('xinzang');
    	$part1_list.=read_CS('qita');
	
	$part1_list.="\\end{longtable} \n注：以上检测疾病中无“*”标识的疾病与基因间的相关性均来自于国内外指南或专家共识；带“*”标识的疾病与基因间的相关性仅来自于当前最新的科学研究，尚未在指南和共识中明确指出。\n\\noindent \\\\ \\\\ \\\\ \n";
		
	}
	
	if($types eq "WES")
	{
		$part1_list.="{\\sym{附表$fb_num：}}";
		$part1_list.=<<'END';
本次检测疾病相关基因如下：\\
\renewcommand\arraystretch{1.8}
\fontsize{10.5}{12.6}\selectfont
\begin{longtable}{ L{7cm}  L{8.2cm} }
\hline
\makecell[c]{\sym{检测疾病}} & \makecell[c]{\sym{相关基因}} \\ \hline
END

		open IN, "$genelist" or die $!;
		while(<IN>)
		{
			chomp;
			@cut=split /\t/,$_;
			if($cut[0] eq $samples )
			{
				$part1_list.="\\makecell[{{p{6.9cm}}}]{$cut[1]} & \\makecell[{{p{8.1cm}}}]{$cut[2]} \\\\ \\hline \n";
			}
		}
		close IN;
		
		$part1_list.=<<'END';
\end{longtable}
\noindent \\
END
	}
}


################### qc #################

my $qc_flag=0;
my ($qc_qc30,$qc_10,$qc_20,$qc_dep);
open IN, $qc_stat or die $!;
while(<IN>)
{
	chomp;
	next if(/^$/);
	@cut=split /\t/,$_;
	if($cut[1] eq "Raw reads")
	{
		$qc_flag=1;

	}
	elsif($cut[1] eq "LibID:")
	{
		$qc_flag=2;
	}
	if($qc_flag==1 and $cut[0] eq $samples){ $qc_qc30="$cut[10];$cut[11]"; $qc_qc30=~s/\%/\\\%/g; }
	if($qc_flag==2 and $cut[0] eq $samples){ $qc_10=$cut[33];$qc_20=$cut[32];$qc_dep=$cut[3]; }
		
}
close IN;

$fb_num++;
$part6="\\fontsize{10.5}{12.6}\\selectfont\n{\\sym{附表$fb_num：}} ";
$part6.=<<'END';
本次检测的主要测序参数如下：\\
本次检测采用方法为目标序列捕获高通量测序技术，采用 Illumina 公司 HiSeq 系列测序平台，产出数据的质量及准确性皆为业内最高标准，配合严谨的信息分析流程和数据解读方法，确保您宝贵的基因信息得到最科学的管理。\\
\vspace*{-3mm}
\renewcommand\arraystretch{1.5}
\begin{center}
\begin{tabular}{  C{3.5cm} | C{3.5cm} | C{3.5cm} | C{3.5cm}  } \hline
\multicolumn{4}{c}{\sym{目标区域高通量测序参数}} \\ \hline
{\sye{Q30}} & {\sye \sym{10X覆盖度（\%）}} & {\sye \sym{20X覆盖度（\%）}} & {\sym{平均测序深度}} \\ \hline
END

$part6.=qq[$qc_qc30 & $qc_10 & $qc_20 & $qc_dep X \\\\ \\hline \n];
$part6.=<<'END';
\end{tabular}
\end{center}
END

################### check #################

my $now = `date "+%F %T"`;
my ($t1,$t2,$t3,$t4,$t5)=$now=~/(\d\d\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):\d+/;
my $times_now="$t1-$t2-$t3";
$check=qq[\\zihao{5}
\\begin{center}
\\begin{tabular}{  C{4.9cm}  C{4.9cm}  C{5cm}  } \\hline
报告分析：$person & 报告审核：冷雪 & 报告日期：$times_now \\\\ \\hline
\\end{tabular}
\\end{center}];



#######################################################
sub read_CS{
	my $csid=shift;
	my $line='';
	my %head_Dic=('xinlv'=>'遗传性心律失常（59个）',
	  'xinji'=>'遗传性心肌病（101个）',
	  'zhudongmai'=>'遗传性主动脉病（12个）',
	  'xinzang'=>'先天性心脏病（24个）',
	  'qita'=>'其他心脏相关疾病（12个）'
	);
	$line="\\multicolumn{2}{c}{\\sym \\sye{$head_Dic{$csid}}} \\\\ \\hline\n\\makecell[c]{\\sym{检测疾病}} & \\makecell[c]{\\sym{相关基因}} \\\\  \\hline\n";
	open IN, "/PUBLIC/pipline/script/Report/XKFW/V2_171130/database/cs_disease_".$csid.".txt" or die $!;
	while(<IN>)
	{
		chomp;
		@cut=split /\t/,$_;
		$line.="\\raisebox{-0.2ex}{\\shortstack[l]{$cut[0]}} & $cut[1] \\\\ \\hline \n";
	}
	close IN;
	
	return $line;
}

sub judge_lastest_sitedb{
        my $product_site=shift;
        my @all_date;
        my $lastest_date;
        my $date;

        my @all_path=glob "$product_site/*";
        foreach(@all_path)
        {
	        $date=(split /\//,$_)[-1];
                if($date=~/2018心康送样信息单-(\d+)\.xlsx/)
                {
                        push @all_date,$1;
                }
        }

        $lastest_date= max(@all_date);
        my $lastest_sitedb_path ="$product_site/2018心康送样信息单-$lastest_date.xlsx";
	my $lastest_sitedb_path2="$product_site/2018心康送样信息单-$lastest_date.txt";
        print "$lastest_sitedb_path\n";
        return $lastest_sitedb_path,$lastest_sitedb_path2;
}


############################# report latex ##############################

$yemei=<<'END';
%-*- coding: UTF-8 -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[UTF8,noindent,twoside]{ctexrep}
\usepackage{blkarray,fancyvrb,ltxtable,multirow,colortbl,array,makecell,fancyhdr,graphicx,multicol,balance,amssymb,setspace,dcolumn,bo
oktabs,tabularx,longtable,wallpaper,tabu,makeidx,algorithm,algorithmic,changepage,fontspec,wrapfig,picinpar,enumerate,pdfpages,textcom
p,supertabular,paralist}
\usepackage{flushend}
\usepackage[abs]{overpic}
\setlength{\parindent}{0pt}
\usepackage[a4paper,left=2cm,right=2cm,top=3cm,bottom=3cm]{geometry}
\usepackage[table]{xcolor}
\usepackage[center]{titlesec}
\graphicspath{{/PUBLIC/pipline/script/Report/XKFW/V2_171130/photo/}}
\definecolor{green}{HTML}{539B34}
\definecolor{yellow}{HTML}{FABE00}
\definecolor{orange}{HTML}{ED6D2A}
\definecolor{red}{HTML}{E60012}
\definecolor{black}{HTML}{000000}
\definecolor{MyDarkBlue}{HTML}{539B35}
\definecolor{putonggreen}{HTML}{44CC6A}
\definecolor{MyLightGreen}{HTML}{F4F7EF}
\definecolor{MyGreen}{HTML}{539B35}
\definecolor{tabcolor}{HTML}{539B35}
\definecolor{MyFontGray}{HTML}{4C4949}
\definecolor{DarkBlue}{HTML}{004271}
\setmainfont{SourceHanSansCN-Light.ttf}
\setsansfont{SourceHanSansCN-Light.ttf}
\setmonofont{SourceHanSansCN-Light.ttf}
\setCJKmainfont{SourceHanSansCN-Light.ttf}
\setCJKsansfont{SourceHanSansCN-Light.ttf}
\setCJKmonofont{SourceHanSansCN-Light.ttf}
\newfontfamily\sye{SourceHanSansCN-Medium.ttf}
\newCJKfontfamily\sym{SourceHanSansCN-Medium.ttf}
%table format command
\newcolumntype{C}[1]{>{\centering\arraybackslash}m{#1}}
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\arraybackslash}m{#1}}
\newcolumntype{Y}{>{\centering\arraybackslash}X}
\newcommand{\tabincell}[2]{\begin{tabular}{@{}#1@{}}#2\end{tabular}}
\titleformat{\section}{\centering\zihao{-2}}{}{2em}{}
\newcommand\mgape[1]{\gape{$\vcenter{\hbox{#1}}$}}  %%%表格中图片居中
%pageHeader setting
\pagestyle{fancy}
%\pagestyle{myheadings}
\fancyhf{}
\renewcommand\headrulewidth{0.6pt}
\renewcommand{\footrulewidth}{0.2pt}
%\setlength{\headheight}{4cm}
END
$yemei.=qq[$header];
$yemei.=qq[$footer];
$yemei.=<<'END';
%\fancyfoot[CO,CE]{\textbf{\thepage}}
\makeatletter %双线页眉
\def\headrule{{\color{black}\if@fancyplain\let\headrulewidth\plainheadrulewidth\fi
\hrule\@height 0.0pt \@width\headwidth\vskip1pt%上面线为1pt粗
\hrule\@height 1.0pt\@width\headwidth %下面0.5pt粗
\vskip-2\headrulewidth\vskip-3pt} %两条线的距离1pt
\vspace{3mm}} %双线与下面正文之间的垂直间距
\begin{document}
%%%%%%%%%%%%%%%%%%%%%%% Report Face %%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage
\vspace*{-4mm}
END

if($types eq "DM")
{
	$biaoti="单基因糖尿病基因检测报告";
}
elsif($types eq "CS")
{
	$biaoti="遗传性心血管病基因检测报告";
}
elsif($types eq "GXY")
{
	$biaoti="单基因高血压/血钾异常基因检测报告";
}
elsif($types eq "WES")
{
        $biaoti="单基因遗传病基因检测报告";
}
$biaoti2=qq[\\section*{\\sym{$biaoti}}
\\zihao{-4} {\\sym{检测编号：$samples}} \\\\ \\\\
\\color{MyFontGray} \\fontsize{10.5}{12.6}\\selectfont
\\renewcommand\\arraystretch{1.5}
\\noindent
\\begin{tabular}{  L{3.5cm}  L{3.5cm}  L{3.5cm}  L{4cm}  } \\hline
];

if($types eq "DM")
{
	$part7_1="1、本次检测重点关注疾病数据库收录的与单基因糖尿病相关的38个基因的全部外显子区域";
}
if($types eq "CS")
{
	$part7_1="1、本次检测重点关注疾病数据库收录的与遗传性心血管病相关的179个基因的全部外显子区域";
}
if($types eq "GXY")
{
	$part7_1="1、本次检测重点关注疾病数据库收录的与单基因高血压/血钾异常病相关的41个基因的全部外显子区域";
}
if($types eq "WES")
{
	$part7_1="1、本次检测重点关注疾病数据库收录的与".$disease_name."疾病相关的".$disease_gene."个基因的全部外显子区域";
}


$part7=<<'END';
\begin{spacing}{1.4}
{\sym{备注：}}\\
END
$part7.=$part7_1."、外显子与内含子衔接处的片段，不包括基因调控区及深度内含子区，检测仅能提供上述覆盖区域内的位点变异和小片段插入/缺失（<25bp）等信息分析。此外，由于二代测序的局限性，现有的技术难以捕获部分GC含量高的基因区域，GC含量高的区域也会有未能覆盖的情况。\\\\ \n";
$part7.=<<'END';
2、本检测不适用于分析基因组结构变异（如大片段拷贝数缺失、插入或倒位）、动态变异、复杂重组等特殊类型的变异。\\ \\
\end{spacing}
END

my $all_tex;
if($negative == 0)
{
	$all_tex=$yemei.$biaoti2.$info_sample.$part1.$part2."\n \\end{longtable}\n\\fontsize{10.5}{12.6}\\selectfont\n".$part3."\\fontsize{10.5}{12.6}\\selectfont \n\\noindent \\\\ \n\\begin{spacing}{1.3}\n{\\sym {检测结果说明：}} \\\\ \n".$part4."\\end{spacing}\n".$part5.$part6_0.$part1_list.$part6.$part7.$part8.$check."\n\\end{document} \n";
}
elsif($negative == 1)
{
	$all_tex=$yemei.$biaoti2.$info_sample.$part1.$part6_0.$part1_list.$part6.$part7.$part8.$check."\n\\end{document} \n";
}
elsif($negative == 2)
{
	$all_tex=$yemei.$biaoti2.$info_sample.$part1.$part1_list.$part6.$part7.$part8.$check."\n\\end{document} \n";
}

open OUT, ">$out/$samples.$types.tex";
print OUT $all_tex;
close OUT;
