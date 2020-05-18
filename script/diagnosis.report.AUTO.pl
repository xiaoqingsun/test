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


my ($samples,$site_file,$out,$types,$qc_stat,$hospital,$usage,$genelist,$feature_file,$cnv_file)=('','','','','','','','','','');
GetOptions(
    'sample=s'=> \$samples,
    'hospital=s'=> \$hospital,
    'site=s' => \$site_file,
    'out=s'  => \$out,
    'type=s' => \$types,
    'qc=s'   => \$qc_stat,
    'genelist=s'=> \$genelist,  #sxq 5.29
    'feature=s' => \$feature_file,
    'cnv=s'  => \$cnv_file,
    #'gene_dis=s' => \$gene_disease_file,
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
my ($name,$sex,$year,$id_len,$shizi,$clinical_info,$jiazu,$yiyuan,$header,$footer)=('','','',0,"外周血","无","无","NULL",'','');
my ($yemei,$biaoti,$biaoti2,$info_sample,$info1,$part1,$part2,$part3,$part4,$part5,$part6,$part6_0,$part7_1,$part7,$part8,$check,$part9,$part10);
my ($part4_1,$part4_2,$part4_3);

my $feature="无临床信息";
my $part1_list='';
my $negative=0;
my $negative2=0;
my $fb_num=0;
#my $gene_disease_file="/PUBLIC/work/panqi/auto_ACMG/OMIM_phen_gene.gene2omim";
#my $bin = "/PUBLIC/work/panqi/diagnosis_report/script";
my $bin = "/PUBLIC/pipline/script/Report/XKFW/V2_171130/script";
my %omim_name;
my %omim_name2;

################### part 0 sample info #####################
my $hutongbiao="/PUBLIC/pipline/database/sheet_hutong/";
my ($file,$file2) = judge_lastest_sitedb($hutongbiao);
my $disease_script= judge_lastest_sitedb2('Disease');
`perl $bin/spectrum_XLSX.pl $file` unless(-e $file2);

=head
open IN, $disease_script or die $!;
while(<IN>)
{
	chomp;
	@cut=split /\t/,$_;
	my $omimid;
	my $omimnm;
	
	next if(/#NewID/);
	if($cut[0]=~/OMIM:(\d+)/)
	{
		$omimid=$1;	
	}
	$omimnm=(split /【/, $cut[1])[0];
	$omim_name{$omimid}=$omimnm;
	$omim_name2{$omimid}=$cut[2];
}
close IN;
=cut

open IN, $feature_file or die $!;
while(<IN>)
{
	chomp;
	@cut=split /\t/, $_;
	if($cut[0] eq $samples && $cut[1])
	{
		$feature=$cut[1];
		#$feature=~s/\\/\\\\/g;
		#$feature=~s/\%/\\\%/g;
	}
}
close IN;

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
		if($year=~/岁/ or $year=~/月/ or $year=~/不详/){print $year;}
		else{$year=$year."岁"}#wph 180823
		$yiyuan=$cut[10];
		$shizi=$cut[11];
		$cut[1]=~s/\//-/g;
		
		my $id_shizi=(split /\//,$cut[3])[0];
		if($id_shizi=~/([0-9]+)/)
		{
			if(length($1) == 8 && $shizi=~/0|None/ )
			{
				$shizi="口腔拭子";
			}
		}
		else{
			if(length($1) != 8 && $shizi=~/0|None/ ){ $shizi="外周血"; }
		}
		#$id_len=length $cut[3];
		#if($id_len > 8){$shizi="外周血";}
		$info1 = $shizi;
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
\\multicolumn{4}{l}{ \\fontsize{9.5}{12} \\selectfont{ \\makecell[{{p{16.5cm}}}]{  {\\sym{临床信息:}}\\\\ $feature }} } \\\\ \\hline
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
	if($yiyuan=~/糖友/ &&  $yiyuan eq $cut[0])
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



################### part 1 CNV info #################

my %het_hom = ('het'=>'杂合','hom'=>'纯合','unknown'=>'未明','hem'=>'半合子');
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
    'splicing SNV'=>'剪切区变异',
    'splicing INDEL' =>'剪切区变异');

#   'splicing SNV'=>'剪切区单核苷酸变异');

my %index;
my (%gene_site,%gene_site_level,%gene_site_level_fb,%disease_site_level,%gene_level,%gene_level_fb,%gene_site_level2,%gene_level2,%gene_level2_fb,%gene_level3,%gene_fb, %disease_fb,%gene_cnv,%gene_cnv2);
my ($disease_name,$disease_gene,$disease_gene_line);
my %disease_gene_line_hash;
my (%site_dis,%site_dis_bz,%site_dis_fb);
my ($site_info,$exon);
my ($c1,$c2,$c3,$c4,$c5,$c6,$c6_2,$c7,$c8,$c8_2,$c9);
my ($sn,$s1,$s2,$s3)=(0,'','','');
my (@paper,@arr,@arr2);
my $num_line=0;
my $num_line_yang=0;
my $num_line_yang2=0;
my $person="NULL";
my %second_line=();
my %site_info_tmp  =();
my %site_info_tmp2 =();
my $level_site;
my %gene_omim;

#open IN, $gene_disease_file or die $!;
open IN, $disease_script or die $!;
while(<IN>)
{
	chomp;
	@cut=split /\t/,$_;
	next if($cut[4] eq "-");
	my $omimnm=(split /【/, $cut[4])[0];
	my $omimtp;

=head
	if(!exists $omim_name{$cut[1]})
	{
		$omimnm='-';
		$omimtp='-';
	}
	else{
		$omimnm=$omim_name{$cut[1]};
		$omimtp=$omim_name2{$cut[1]};
	}

	if($cut[0] eq "CFHR3")
	{
		$cut[0]="CFHR3";
	}
	
#	$gene_cnv{$cut[0]} ="$omimnm($cut[4])[$omimtp]"     if(!exists $gene_cnv{$cut[0]});
#	$gene_cnv{$cut[0]}.="\\\\$omimnm($cut[4])[$omimtp]" if( exists $gene_cnv{$cut[0]});
=cut

	if(!exists $gene_omim{$cut[3]}{$cut[0]}){
		$gene_omim{$cut[3]}{$cut[0]} = 1;
		if($cut[7] eq "-" && $cut[4] ne "-"){$cut[7]="不明";}
		if(!exists $gene_cnv{$cut[3]}){
			$gene_cnv{$cut[3]} ="$omimnm($cut[7])";
		}else{
			$gene_cnv{$cut[3]}.="\\\\$omimnm($cut[7])";
		}
	}
}
close IN;


################### part 2 site info #################
#open IN, "<:utf8", "$site_file" or die $!;

open IN, $site_file or die $!;
while(<IN>)
{
	chomp;
	$_=~s/\_/\\\_/g;
	$_=~s/\%/\\\%/g;
	@cut=split /\t/,$_;
	$level_site=0;
	my ($gene_info, $dis_name, $dis_info, $acmg_info, $site_miaoshu_info);
	
	if($num_line == 0)
	{
		for(my $i=0;$i<@cut;$i++)
		{
			$index{$cut[$i]}=$i;
		}		
	}
	else{
		next unless ($cut[0]=~/$samples/);
		next if($cut[0]=~/^NN/);
		if($cut[0]=~/^####/){$person=$cut[$index{'解读人'}] if(exists $index{'解读人'} && $cut[$index{'解读人'}] ne "-");next;}
		#$person=$cut[$index{'解读人'}] if(exists $index{'解读人'} && ($cut[$index{'解读人'}] ne "-") && ($cut[$index{'解读人'}] ne "") && $cut[$index{'解读人'}]);
		$person="NULL";  #if(!exists $index{'解读人'});
		$site_info=join("\_",@cut[1,2,4,5]);

		if($site_info eq "-_-_-_-")
		{
			$negative=2;
			next;
		}

		if($cut[$index{'PDF\_NMID'}] eq "\.\/\/")
		{
			$exon="NULL";
		}else{
			$exon=join "\\\\", ( split /:/, $cut[$index{'PDF\_NMID'}] )[0,1];
		}
		
		if($cut[$index{'基因描述'}]=~/\[auto\]/ || $cut[$index{'基因描述'}]=~/\[done\]/ || $cut[$index{'基因描述'}]=~/\[check\]/ )
		{
			$gene_info=(split /\[auto\]/, $cut[$index{'基因描述'}])[0] if($cut[$index{'基因描述'}]=~/\[auto\]/);
			$gene_info=(split /\[done\]/, $cut[$index{'基因描述'}])[0] if($cut[$index{'基因描述'}]=~/\[done\]/);
			$gene_info=(split /\[check\]/, $cut[$index{'基因描述'}])[0] if($cut[$index{'基因描述'}]=~/\[check\]/);
		}
		else{	$gene_info= $cut[$index{'基因描述'}]; }
		$gene_info=~s/gt;/\>/g;
		$gene_info=~s/lt;/\</g;
		
		if($cut[$index{'疾病名称'}]=~/\[auto\]/ || $cut[$index{'疾病名称'}]=~/\[done\]/ || $cut[$index{'疾病名称'}]=~/\[check\]/ || $cut[$index{'疾病名称'}]=~/\[product\]/)
		{
			$dis_name=(split /\[auto\]/, $cut[$index{'疾病名称'}])[0] if($cut[$index{'疾病名称'}]=~/\[auto\]/);
			$dis_name=(split /\[done\]/, $cut[$index{'疾病名称'}])[0] if($cut[$index{'疾病名称'}]=~/\[done\]/);
			$dis_name=(split /\[check\]/, $cut[$index{'疾病名称'}])[0] if($cut[$index{'疾病名称'}]=~/\[check\]/);
			$dis_name=(split /\[product\]/, $cut[$index{'疾病名称'}])[0] if($cut[$index{'疾病名称'}]=~/\[product\]/);
		}
		else{	$dis_name= $cut[$index{'疾病名称'}]; }
		$dis_name=~s/gt;/\>/g;
		$dis_name=~s/lt;/\</g;

		if($cut[$index{'疾病描述'}]=~/\[auto\]/ || $cut[$index{'疾病描述'}]=~/\[done\]/ || $cut[$index{'疾病描述'}]=~/\[check\]/ || $cut[$index{'疾病描述'}]=~/\[product\]/ )
		{
			$dis_info=(split /\[auto\]/, $cut[$index{'疾病描述'}])[0] if($cut[$index{'疾病描述'}]=~/\[auto\]/);
			$dis_info=(split /\[done\]/, $cut[$index{'疾病描述'}])[0] if($cut[$index{'疾病描述'}]=~/\[done\]/);
			$dis_info=(split /\[check\]/, $cut[$index{'疾病描述'}])[0] if($cut[$index{'疾病描述'}]=~/\[check\]/);
			$dis_info=(split /\[product\]/, $cut[$index{'疾病描述'}])[0] if($cut[$index{'疾病描述'}]=~/\[product\]/);
		}
		else{	$dis_info= $cut[$index{'疾病描述'}]; }
		$dis_info=~s/gt;/\>/g;
		$dis_info=~s/lt;/\</g;

		if( $cut[$index{'位点描述'}]=~/\[auto\]/ || $cut[$index{'位点描述'}]=~/\[done\]/ || $cut[$index{'位点描述'}]=~/\[check\]/ )
		{
			$site_miaoshu_info=(split /\[auto\]/, $cut[$index{'位点描述'}])[0] if($cut[$index{'位点描述'}]=~/\[auto\]/);
			$site_miaoshu_info=(split /\[done\]/, $cut[$index{'位点描述'}])[0] if($cut[$index{'位点描述'}]=~/\[done\]/);
			$site_miaoshu_info=(split /\[check\]/, $cut[$index{'位点描述'}])[0] if($cut[$index{'位点描述'}]=~/\[check\]/);
		}
		else{	$site_miaoshu_info= $cut[$index{'位点描述'}]; }
		$site_miaoshu_info=~s/gt;/\>/g;
		$site_miaoshu_info=~s/lt;/\</g;
		$cut[$index{'位点描述'}]= $site_miaoshu_info;

		my $acmg_line;
		$acmg_line=$cut[$index{'ACMG条目'}]     if(exists $index{'ACMG条目'});
		$acmg_line=$cut[$index{'ACMG对应条目'}] if(exists $index{'ACMG对应条目'});
				    
		my @acmg_array = split /\+/,$acmg_line;
		#@acmg_array = split /\+/,$acmg_line if($acmg_line=~/\+/);
		@acmg_array = split /,/,$acmg_line if($acmg_line=~/,/);
		@acmg_array = split /，/,$acmg_line if($acmg_line=~/，/);
		@acmg_array = split / /,$acmg_line if($acmg_line=~/ /);
		my @acmg_array_sort=sort @acmg_array;
		$acmg_info=$acmg_array_sort[0] if(@acmg_array_sort==1);
		$acmg_info=join(" ", @acmg_array_sort[0..$#acmg_array_sort]) if(@acmg_array_sort>1);
		print STDERR $acmg_info."\n";
		
		if($cut[0]=~/^N0/ or !exists $gene_site{$cut[7]}{$site_info})
		{		
			$gene_site{$cut[7]}{$site_info}="$cut[$index{'HGNC基因名'}]\t$cut[$index{'CDS'}]\t$cut[$index{'PEP'}]\t$exon\t$het_hom{$cut[31]}\t$mute_des{$cut[9]}\t$cut[21]\t$gene_info" if(exists $index{'HGNC基因名'});
			$gene_site{$cut[7]}{$site_info}="$cut[$index{'基因名称'}]\t$cut[$index{'CDS'}]\t$cut[$index{'PEP'}]\t$exon\t$het_hom{$cut[31]}\t$mute_des{$cut[9]}\t$cut[21]\t$gene_info" if(exists $index{'基因名称'});
		}

		if($cut[0]=~/^###/)
		{
			push @{ $site_dis_bz{$cut[7]}{$site_info} }, [ $dis_name, $cut[$index{'遗传模式'}], $acmg_info, $cut[$index{'致病性结论'}] ]  if(exists $index{'ACMG条目'});
			push @{ $site_dis_bz{$cut[7]}{$site_info} }, [ $dis_name, $cut[$index{'遗传模式'}], $acmg_info, $cut[$index{'位点结论'}] ]  if(exists $index{'ACMG对应条目'});
		}
		elsif($cut[0]=~/^##/)
		{
			$negative2=1;
			if(exists $index{'ACMG条目'})
			{
				#$gene_fb{ $cut[$index{'基因名称'}] } = $cut[$index{'基因描述'}];
				$gene_fb{ $cut[$index{'基因名称'}] } = $gene_info;
				$disease_fb{ $cut[$index{'基因名称'}] }{ $dis_name } = "$dis_name : $dis_info";
				if(!exists $gene_site_level_fb{$cut[7]}{$site_info} )
				{
					$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{ $cut[$index{'致病性结论'}] };
				}
				elsif($gene_site_level_fb{$cut[7]}{$site_info} > $level_value{ $cut[$index{'致病性结论'}]} ){
					$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{ $cut[$index{'致病性结论'}] };
				}
				#$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{ $cut[$index{'致病性结论'}] };
				push @{$gene_level_fb{$cut[7]}},$level_value{ $cut[$index{'致病性结论'}] };
				push @{ $site_dis_fb{$cut[7]}{$site_info}{$acmg_info} }, [ $dis_name, $cut[$index{'遗传模式'}], $acmg_info, $cut[$index{'致病性结论'}] ];
			}
			if(exists $index{'ACMG对应条目'})
			{
				#$gene_fb{ $cut[$index{'HGNC基因名'}] } = $cut[$index{'基因描述'}];
				$gene_fb{ $cut[$index{'HGNC基因名'}] } = $gene_info;
				$disease_fb{ $cut[$index{'HGNC基因名'}] }{ $dis_name } = "$dis_name : $dis_info";
				if(!exists $gene_site_level_fb{$cut[7]}{$site_info} )
				{
					$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{ $cut[$index{'位点结论'}] };
				}
				elsif($gene_site_level_fb{$cut[7]}{$site_info} > $level_value{ $cut[$index{'位点结论'}]} ){
					$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{ $cut[$index{'位点结论'}] };
				}
				#$gene_site_level_fb{$cut[7]}{$site_info}=$level_value{ $cut[$index{'位点结论'}] };
				push @{$gene_level_fb{$cut[7]}},$level_value{$cut[$index{'位点结论'}]};
				push @{ $site_dis_fb{$cut[7]}{$site_info}{$acmg_info} }, [ $dis_name, $cut[$index{'遗传模式'}], $acmg_info, $cut[$index{'位点结论'}] ];
			}
			if($cut[$index{'参考文献'}] ne '-' && $cut[$index{'参考文献'}] ne '.')
			{
				push @paper, split(/\|/, $cut[$index{'参考文献'}]);
			}
		}
		else{
			if(exists $index{'位点结论'})
			{
				$level_site=$level_value{ $cut[$index{'位点结论'}]  };
				push @{ $gene_site_level{$cut[7]}{$site_info} }, $level_value{ $cut[$index{'位点结论'}]  };
				if(!exists $disease_site_level{$cut[7]}{$site_info}{$dis_name})
				{
					push @{ $site_dis{$cut[7]}{$site_info}{ $acmg_info } }, [ $cut[$index{'位点结论'}], $dis_name, $cut[$index{'遗传模式'}], $dis_info, $acmg_info, $cut[$index{'位点描述'}] ];
					$disease_site_level{$cut[7]}{$site_info}{$dis_name}=$cut[$index{'位点结论'}];
				}
			}
			elsif(exists $index{'致病性结论'})
			{
				$level_site=$level_value{ $cut[$index{'致病性结论'}]  };
				push @{ $gene_site_level{$cut[7]}{$site_info} }, $level_value{ $cut[$index{'致病性结论'}]  };
				if(!exists $disease_site_level{$cut[7]}{$site_info}{$dis_name})
				{
					push @{ $site_dis{$cut[7]}{$site_info}{ $acmg_info } }, [ $cut[$index{'致病性结论'}], $dis_name, $cut[$index{'遗传模式'}], $dis_info, $acmg_info, $cut[$index{'位点描述'}] ];
					$disease_site_level{$cut[7]}{$site_info}{$dis_name}=$cut[$index{'致病性结论'}];
				}
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
			if(!exists $site_info_tmp2{$site_info} && $level_site<3)
			{
				$num_line_yang2++;
				$site_info_tmp2{$site_info}=1;
			}
			#$person=$cut[$index{'解读人'}] if(exists $index{'解读人'});
			#$person="NULL" if(!exists $index{'解读人'});
		}
	}
	$num_line++;
}
close IN;

unless(%gene_site_level)
{
	$negative=1 if($negative2==1);
	$negative=2 if($negative2==0);
}
if(%gene_site_level)
{
	$negative=3 if($negative2==0);
}

################### splicing and nonsense mutation #################

my (%splice_site, %stopgain_site);
my ($splice_site_num, $stopgain_site_num)=(0,0);
my $dis_special;
my $site_file2 = $site_file;
$site_file2=~s/IDT\.pdf\.txt/IDT\.right\.2\.txt/;

open IN, $site_file2 or die $!;
while(<IN>)
{
	chomp;
	$_=~s/\_/\\\_/g;
	$_=~s/\%/\\\%/g;
	@cut=split /\t/,$_;
	next unless ($cut[1]=~/$samples/);

	$site_info=join("\_",@cut[2,3,5,6]);
	
	next if(exists $gene_site{$cut[8]}{$site_info});
	next if($cut[10]!~/splicing SNV/ and $cut[10]!~/stopgain/);
	
	if($cut[14] eq "\.\/\/" or $cut[14] eq "-")
	{
		$exon="NULL";
	}else{
		$exon=join "\\\\", ( split /:/, $cut[14] )[0,1];
	}
	
	#@arr = split /\|\|/,$cut[84];
	#@arr2= split /\|\|/,$cut[85];
	#for(my $i=0; $i<@arr; $i++)
	#{
	#	$arr2[$i]='-' if($arr2[$i] eq '');
	#	$dis_special ="$arr[$i]($arr2[$i])" if($i==0);
	#	$dis_special.="\\\\$arr[$i]($arr2[$i])" if($i>0 );
	#}

	if(exists $gene_cnv{$cut[8]}){
		$dis_special = $gene_cnv{$cut[8]};
	}
	else{
		$dis_special = "-(-)";
	}
	
	if($cut[10]=~/splicing SNV/)
	{
		if($cut[15]=~/c\.(\d+)([\+\-])(\d+)([ATCGN]+)>([ATCGN]+)/ )
		{
			if($3<3){
				$splice_site{$site_info} = "$cut[8]\t$cut[15]\t$cut[16]\t$exon\t$het_hom{$cut[32]}\t$mute_des{$cut[10]}\t$dis_special"; 
				$splice_site_num++;
			}
		}
	}
	if($cut[10]=~/stopgain/)
	{
		if($cut[5] ne "-" && $cut[6] ne "-" && length($cut[5])==1 && length($cut[6])==1)
		{
			$stopgain_site{$site_info} = "$cut[8]\t$cut[15]\t$cut[16]\t$exon\t$het_hom{$cut[32]}\t$mute_des{$cut[10]}\t$dis_special";
			$stopgain_site_num++;
		}
	}

	$gene_site{$cut[8]}{$site_info}=1;
}
close IN;

if($splice_site_num==0 && $stopgain_site_num==0)
{
	$part9="";
}
elsif($splice_site_num>0 || $stopgain_site_num>0)
{
	$part9=<<'END';
%\noindent \\ \\ \\
\newpage
\begin{spacing}{1.3}
\fontsize{10.5}{12.6}\selectfont
END
	if($negative<=1)
	{
		$part9.=qq[{\\sym{四、特殊类型变异检测结果：}} \\\\ \n];
	}
	elsif($negative>1){
		$part9.=qq[{\\sym{三、特殊类型变异检测结果：}} \\\\ \n];
	}

	$part9.=qq[\\renewcommand\\arraystretch{1.3} \n];
	$part9.=qq[\\fontsize{10.5}{12.6}\\selectfont{\\color{MyFontGray} 受检者检测到 $splice_site_num 个剪切位点±1/2位置的变异、$stopgain_site_num 个蛋白合成提前终止的变异，由于目前研究较少，其致病性无法确定，但是该类型突变可能会影响基因及基因产物功能，请结合临床综合考虑，具体信息如下表所示：} \n \\end{spacing} \n];

	if($splice_site_num>0)
	{
		$part9.=<<'END';
\\ \\
\fontsize{10.5}{12.6}\selectfont{\sym{剪切变异：}}
\fontsize{8}{9.6}\selectfont
%\begin{longtable}{| C{15mm} | C{16mm} | C{18mm} | C{17mm} | C{8mm} | C{8mm} | C{18mm} | L{35mm} | }
%\begin{longtable}{| C{13mm} | C{16mm} | C{17mm} | C{16mm} | C{7mm} | C{8mm} | m{19mm}<{\centering} | m{34mm} | m{10mm} | }
\begin{longtable}{| C{15mm} | C{18mm} | C{19mm} | C{22mm} | C{7mm} | C{7mm} | m{40mm} | }
\hline
{\sym{变异基因}} & {\sym{核苷酸变异}} & {\sym{氨基酸变异}} & \makecell[c]{ \sym{转录本}\\ \sym{外显子} } & {\sym{变异状态}} & {\sym{变异类型}} & \makecell[c]{\sym{相关疾病(遗传模式)}}  \\ \hline
END

		foreach my $key1 (sort {$a cmp $b} keys %splice_site)
		{
			($c1,$c2,$c3,$c4,$c5,$c6,$c7)=split /\t/, $splice_site{$key1};
			my $c11=word_wrap($c1,10);
			$part9.=qq[ $c11 & \\makecell[c]{$c2} & \\makecell[c]{$c3} & \\makecell[c]{$c4} & $c5 & \\makecell[{}{m{5.8mm}}]{$c6} & \\makecell[{}{m{36mm}}]{$c7} \\\\ \\hline \n];
		}
		$part9.="\\end{longtable}\n";
	}

	if($stopgain_site_num>0)
	{		
		$part9.=<<'END';
%\vspace*{3mm}
\\ \\ 
%\renewcommand\arraystretch{1.3}
\fontsize{10.5}{12.6}\selectfont{\sym{无义变异：}}
\fontsize{8}{9.6}\selectfont
%\begin{longtable}{| C{15mm} | C{16mm} | C{18mm} | C{17mm} | C{8mm} | C{8mm} | C{18mm} | L{35mm} | }
%\begin{longtable}{| C{13mm} | C{16mm} | C{17mm} | C{16mm} | C{7mm} | C{8mm} | m{19mm}<{\centering} | m{34mm} | m{10mm} | }
\begin{longtable}{| C{15mm} | C{18mm} | C{19mm} | C{22mm} | C{7mm} | C{7mm} | m{40mm} | }
\hline
{\sym{变异基因}} & {\sym{核苷酸变异}} & {\sym{氨基酸变异}} & \makecell[c]{ \sym{转录本}\\ \sym{外显子} } & {\sym{变异状态}} & {\sym{变异类型}} & \makecell[c]{\sym{相关疾病(遗传模式)}}  \\ \hline
END

		foreach my $key1 (sort {$a cmp $b} keys %stopgain_site)
		{
			($c1,$c2,$c3,$c4,$c5,$c6,$c7)=split /\t/, $stopgain_site{$key1};
			$part9.=qq[ $c1 & \\makecell[c]{$c2} & \\makecell[c]{$c3} & \\makecell[c]{$c4} & $c5 & \\makecell[{}{m{5.8mm}}]{$c6} & \\makecell[{}{m{36mm}}]{$c7} \\\\ \\hline \n];
		}
		$part9.="\\end{longtable}\n";
	}
}

###################  CNV  ##############
my $part10_2="";
my $num_cnv = 0;
my $num_cnv_gene=0;
my ($cnv_loc,$cnv_copy);
my $gene_cnv_dis;

open IN, $cnv_file or die $!;
while(<IN>)
{
	chomp;
	@cut=split /\t/,$_;
	next if($cut[0] eq "");   ### panqi 190816
	#next if(!exists $gene_cnv{$cut[5]});
	$gene_cnv2{$cut[14]} = 1;
	$num_cnv++;
	
	$cnv_loc="chr".$cut[0].":$cut[1]-$cut[2]";
	
	$cnv_copy=4 if($cut[7]>=1.9 );
	$cnv_copy=3 if($cut[7]>=1.4  && $cut[7]<1.9);
	$cnv_copy=1 if($cut[7]>=0.15 && $cut[7]<0.6);
	$cnv_copy=0 if($cut[7]<0.15 );

	if(exists $gene_cnv{$cut[14]}){
		$gene_cnv_dis=$gene_cnv{$cut[14]}; 
	}
	else{
		$gene_cnv_dis="-(-)"; 
	}
	$part10_2.=qq[ $cnv_loc & $cnv_copy & $cut[14] & \\makecell[{}{m{74mm}}]{$gene_cnv_dis} \\\\ \\hline \n];
}
close IN;

$num_cnv_gene=keys %gene_cnv2;

if($num_cnv==0)
{
	$part10="";
}
else{
	$part10=<<'END';
%\noindent \\ \\ \\
\newpage
\fontsize{10.5}{12.6}\selectfont
END

	if( $negative<=1 && $part9 ne "")
	{
		$part10.=qq[{\\sym{五、染色体拷贝数异常检测结果：}} \\\\ \n];
	}
	elsif($negative<=1 && $part9 eq "")
	{
		$part10.=qq[{\\sym{四、染色体拷贝数异常检测结果：}} \\\\ \n];
	}
	elsif($negative>1 && $part9 ne "")
	{
		$part10.=qq[{\\sym{四、染色体拷贝数异常检测结果：}} \\\\ \n];
	}
	elsif($negative>1 && $part9 eq "")
	{
		$part10.=qq[{\\sym{三、特殊类型变异检测结果：}} \\\\ \n];
	}

	$part10.=qq[\\renewcommand\\arraystretch{1.3} \n];
	$part10.=qq[\\fontsize{10.5}{12.6}\\selectfont{\\color{MyFontGray} 该样本共检测到 $num_cnv 个拷贝数变异，包括 $num_cnv_gene 个基因，详见下表：}\\\\];
	$part10.=<<'END';
\fontsize{9}{10.8}\selectfont
%\begin{longtable}{| C{15mm} | C{16mm} | C{18mm} | C{17mm} | C{8mm} | C{8mm} | C{18mm} | L{35mm} | }
%\begin{longtable}{| C{13mm} | C{16mm} | C{17mm} | C{16mm} | C{7mm} | C{8mm} | m{19mm}<{\centering} | m{34mm} | m{10mm} | }
\begin{longtable}{| C{44mm} | C{10mm} | C{12mm} | m{75mm} | }
\hline
{\sym{染色体位置(hg19)}} & {\sym{拷贝数}} & {\sym{该片段内包含基因}} & \makecell[c]{\sym{相关疾病}\\ \sym{(遗传模式)} } \\ \hline
END

	$part10.=$part10_2."\\end{longtable}\n";

	$part10.=<<'END';
\vspace*{2mm}
\begin{spacing}{1.3}
\zihao{6} \color{MyFontGray}
{\sym{备注：}} \\
1. 拷贝数列{\sye{`1'}}意味着该片段在该患者样本中存在1个拷贝，可能存在片段缺失情况，{\sye{`3'}}意味着该片段在患者样本中存在3个拷贝，可能存在片段重复的情况，{\sye{`0'}}意味着2个拷贝全部缺失。\\
2. 通过适用于外显子拷贝数分析软件（cn.MOPS）的软件对该样本外显子测序结果进行拷贝数分析，由于方法的局限性，检测结果仅供参考，如果该患者需要对拷贝数异常进行检测，建议使MLPA等方法进行进一步检测，以便保证结果的准确性。\\
\end{spacing}
%\vspace*{3mm}
%\begin{spacing}{1.4}
%\zihao{-4}{\sym {\underline {进一步检测建议：} } }\\ \\
%\fontsize{9}{10.8}\selectfont
%如果送检医生认为本报告内容不足以解释临床表型，可考虑进行进一步检测：\\
%1. PANEL 检测：全外显子测序在检测覆盖度方面和定制的 PANEL 检测有一定差距，可以考虑针对某类疾病进行更为细致的 PANEL 检测。\\
%2. 线粒体测序：全外显子检测无法检测到线粒体基因变异。如考虑到线粒体基因变异致病的可能性，可进一步进行线粒体基因检测。\\
%3. CNV 检测：全外显子检测无法非常准确、特异地检测大片段变异（ CNV）。如考虑大片段变异的可能性，可进一步进行 CNV 相关检测。\\
%4. 表观遗传学检测：外显子变异无法展示由于甲基化、乙酰化等表观遗传学改变对基因及基因产物的影响，可进一步进行表观遗传学检测。\\
%\end{spacing}
END

}

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
\vspace*{5mm}
\fontsize{10.5}{12.6}\selectfont
{{\sym\sye 附表 3 ： }}参考文献 \\
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
$part3="\\fontsize{9}{10.8}\\selectfont\n{\\sym{备注：}} ";
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
		if($arr2[$i][1] eq "-"){$arr2[$i][1]="不明";}
		if(!exists $arr_tmp{$key1}{$arr2[$i][0]})
		{
			$part3.="$arr2[$i][0]($arr2[$i][1])" if($i == 0);
			$part3.="、$arr2[$i][0]($arr2[$i][1])" if($i > 0);
			$arr_tmp{$key1}{$arr2[$i][0]}=1;
		}
	}
	$part3.="疾病相关，但目前研究较少，如已经出现该病的临床表型，请结合临床考虑该位点的致病性。";
        @arr2=();
}
$part3.="\\\\ \n";
}

$part3.=<<'END';
\begin{spacing}{1.3}
\zihao{6} \color{MyFontGray}
{\sym{备注：}} \\
{\sym{1、致病性等级}} \\
\hspace*{3mm} 1） {\sym{致病：}} 根据 ACMG 基因变异解读标准，该位点突变有多条致病类证据，是致病的重要原因。\\
\hspace*{3mm} 2） {\sym{可能致病：}} 根据 ACMG 基因变异解读标准，该位点突变有一定数量的致病类证据，其导致疾病的可能性大于 90\%。\\
\hspace*{3mm} 3） {\sym{临床意义未明 1 级：}} 根据 ACMG 基因变异解读标准，该位点突变有一条非常强的致病类证据，但总体证据条目不足，因此其致病性无法准确判定，需结合临床信息综合判断。\\
\hspace*{3mm} 4） {\sym{临床意义未明 2 级：}} 根据 ACMG 基因变异解读标准，该位点突变有一条中等的致病类证据，但总体证据条目不足，因此其致病性无法准确判定，需结合临床信息综合判断。\\
\hspace*{3mm} 5） {\sym{临床意义未明 3 级：}} 根据 ACMG 基因变异解读标准，该位点突变有一条较弱的致病类证据，但总体证据条目不足，因此其致病性无法准确判定，需结合临床信息综合判断。\\
\hspace*{3mm} 6） {\sym{临床意义未明 4 级：}} 根据 ACMG 基因变异解读标准，该位点突变找不到任何证据条目，或证据存在矛盾情况，因此其致病性无法准确判定，需结合临床信息综合判断。\\
\hspace*{3mm} 7） {\sym{临床意义未明 5 级：}} 根据 ACMG 基因变异解读标准，该位点突变有一条较弱的良性类证据，但总体证据条目不足，因此其致病性无法准确判定，需结合临床信息综合判断。\\
\hspace*{3mm} 8） {\sym{可能良性：}} 根据 ACMG 基因变异解读标准，该位点突变有一定数量的良性类证据，其不会致病的可能性大于 90\%。\\
\hspace*{3mm} 9） {\sym{良性：}} 根据 ACMG 基因变异解读标准，该位点突变有多条良性类证据，因此为良性变异，不会导致疾病发生。\\
{\sym\sye{2、ACMG证据条目}} \\
\hspace*{3mm} 1） {\sye{PS:}} 强致病证据     \\
\hspace*{3mm} 2） {\sye{PM:}} 中等致病证据   \\
\hspace*{3mm} 3） {\sye{PP:}} 支持性致病证据 \\
\end{spacing}
\\ \\
END

$part2="";
if(%gene_site_level)
{
	$part1=<<'END';
\noindent \\ \\ \\ 
\begin{spacing}{1.2}
{\sym{二、主要变异检测结果：}} \\
\renewcommand\arraystretch{1.5}
END
	$part1.=qq[\\fontsize{10.5}{12.6}\\selectfont{\\color{MyFontGray} 受检者检测到 $num_line_yang 个与临床表现相关的变异位点，根据ACMG指南进行该位点的综合分析，发现其中 $num_line_yang2 个位点具有致病的可能性，请结合临床综合考虑，具体信息如下表所示：}\\\\];
	$part1.=<<'END';
\end{spacing}
\fontsize{8}{9.6}\selectfont
%\begin{longtable}{| C{15mm} | C{16mm} | C{18mm} | C{17mm} | C{8mm} | C{8mm} | C{18mm} | L{35mm} | }
\begin{longtable}{| C{13mm} | C{16mm} | C{17mm} | C{17mm} | C{7mm} | C{7mm} | m{10mm}<{\centering} | m{18mm}<{\centering} | m{34mm} | }
\hline
{\sym{变异基因}} & {\sym{核苷酸变异}} & {\sym{氨基酸变异}} & \makecell[c]{ \sym{转录本}\\ \sym{外显子} } & {\sym{变异状态}} & {\sym{变异类型}} & \makecell[c]{\sye{ACMG} \\ \sym{条目}}  & {\sym{致病性}} & \makecell[c]{\sym{相关疾病(遗传模式)}} \\ \hline
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
			@arr2=split /\_/,$key2;
			($c1,$c2,$c3,$c4,$c5,$c6,$c6_2,$s2)=split /\t/, $gene_site{$key1}{$key2};
			my $c22=word_wrap($c2,12);
			my $c33=word_wrap($c3,12);
			#$gene_level3{$key1}{g}="{\\sym{基因描述：}} ".$s2."\\\\ \\\\ \n{\\sym{位点描述：}} \\\\ ";
			$gene_level3{$key1}{g}="$key1: $s2 \\\\ \\begin{spacing}{0.8} \\end{spacing} \n";
			$s3='';
			my ($s3_1,$s3_2)=('','');
			
			$sn++;
			my $j=0;
			my $part2_2='';
			my $part2_flag=0;
			foreach my $key3 (sort { $site_dis{$key1}{$key2}->{$a} <=> $site_dis{$key1}{$key2}->{$b} } keys %{ $site_dis{$key1}{$key2} })
			{
				$part2_flag=0;
				my @arr=@{ $site_dis{$key1}{$key2}{$key3} };
				$c8='';
				$c9='';

				for(my $i=0;$i<@arr;$i++)
				{
					#next if($arr[$i][1] ne $key3);
					if($i>0 && $arr[$i][0] ne $arr[0][0]){ $part2_flag=1; }
					$c7   = $arr[$i][0];

					if($arr[$i][2] eq "-" && $arr[$i][1] ne "-")
					{
						$arr[$i][2]="不明";
					}
					$c8="$arr[$i][1]($arr[$i][2])" if($i==0);
					$c8.="\\\\$arr[$i][1]($arr[$i][2])"             if($i>0 && $part2_flag==0);
					$c8.="\\\\$arr[$i][1]($arr[$i][2],$arr[$i][0])" if($i>0 && $part2_flag==1);
					$c9=$key3; ###$arr[$i][4];			
		
					$arr[$i][3]=~s/\~/--/g if($arr[$i][3]=~/\~/);
					if($j==0 && $i==0)
					{
						$gene_level3{$key1}{d}{$arr[$i][1]}=qq[$arr[$i][1]：$arr[$i][3] \\\\ \\begin{spacing}{0.8} \\end{spacing} \n] if(!exists $gene_level3{$key1}{d}{$arr[$i][1]});
#						$s1=
#qq[{\\sym {位点$sn：}} \\\\
#{\\sym{疾病描述：}} \\\\
#$arr[$i][1]：$arr[$i][3] \\\\ \\\\ \n];
					}
					else
					{
						$gene_level3{$key1}{d}{$arr[$i][1]}=qq[$arr[$i][1]：$arr[$i][3] \\\\ \\begin{spacing}{0.8} \\end{spacing} \n] if(!exists $gene_level3{$key1}{d}{$arr[$i][1]});		
					#	$s1.=qq[$arr[$i][1]：$arr[$i][3] \\\\ \\\\ \n];
					}
					
					if($i==0)
					{
						$s3_1="$arr[$i][1]";
						$s3_2=qq[$arr[$i][5] \\\\ \\begin{spacing}{0.8} \\end{spacing} \n];
					}
					else{
						$s3_1.=", $arr[$i][1]";
					}
					#$part2_2.=qq[ $c1 & \\makecell[c]{$c2} & \\makecell[c]{$c3} & \\makecell[c]{$c4} & $c5 & $c6 & $c7 & \\makecell[{}{m{34mm}}]{$c8_2} \\\\ \\hline\n ];
				}
				
				my $max_row=keys %{$site_dis{$key1}{$key2}};
				$s3.="$c1, $c2, $c3 \\\\ \n".$s3_1.": ".$s3_2;
				
				#$part2_2=qq[ & $c7 & \\makecell[{}{m{34mm}}]{$c8} & \\makecell[{}{m{9.5mm}}]{$c9} \\\\ \\cline{7-9}\n ] if($j==0 && $max_row>1);
				
				$part2_2=qq[ & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8}  \\\\ \\cline{7-9}\n ] if($j==0 && $max_row>1);
				$part2_2=qq[ & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8}  \\\\ \\hline\n ] if($j==0 && $max_row==1);
				$part2_2.=qq[ & & & & & & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8} \\\\ \\cline{7-9}\n ] if($j>0 && $j+1!=$max_row);
				$part2_2.=qq[ & & & & & & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8} \\\\ \\hline\n ] if($j>0 && $j+1==$max_row);
				
				$j++;
			}
	#		$part4.=$s1."{\\sym{基因描述：}} ".$s2."\\\\ \\\\ \n{\\sym{位点描述：}} \\\\ ".$s3;
			$gene_level3{$key1}{s}{$key2}=$s3;
			if($j==1){
				$part2.=qq[ $c1 & $c22 & $c33 & \\makecell[c]{$c4} & $c5 & \\makecell[{}{m{5.8mm}}]{$c6} ].$part2_2;
			}else{
				$part2.=qq[ \\multirow{$j}{*}{$c1} & \\multirow{$j}{*}{$c22} & \\multirow{$j}{*}{$c33} & \\multirow{$j}{*}{\\makecell[c]{$c4} } & \\multirow{$j}{*}{$c5} & \\multirow{$j}{*}{ \\makecell[{}{m{5.8mm}}]{$c6} } ].$part2_2;
			}

			$part5.=
qq[{\\sym\\sye{Sanger验证：}} \\\\
$c1:  $c2, $c3 \\\\ ];

			my $fig="$out/fig/$samples-$arr2[0]-$arr2[1].jpg";
			if(-e $fig)
			{
				$part5.=qq[ \\\\ \\centerline{\\includegraphics[width=14cm,height=3.5cm]{$out/fig/$samples-$arr2[0]-$arr2[1].jpg}} \\\\ \\\\ \\\\ \\\\ \n ];
			}
			else{
				$part5.=
qq[参照序列 \\\\
样本序列 \\\\
\\vspace*{4cm}
\\\\];
			}
		}
		
			#$part5="\\\\ \n";  ##### 一代验证 #####

		foreach my $key4 ( keys %{$gene_level3{$key1}{d} } )
		{
			if(!exists $gene_level3{$key4}{disease})
			{
				$part4_1.=$gene_level3{$key1}{d}{$key4};
				$gene_level3{$key4}{disease}=1;
			}
		}
		$part4_2.=$gene_level3{$key1}{g};
		foreach my $key4 ( keys %{$gene_level3{$key1}{s} } )
		{
			$part4_3.=$gene_level3{$key1}{s}{$key4};
		}
	}	
	#$gene_level3{$key1}{g}="{\\sym{基因描述：}} ".$s2."\\\\ \\\\ \n{\\sym{位点描述：}} \\\\ ";
	$part4="{\\sym{疾病描述：}} \\\\ \n".$part4_1."{\\sym{基因描述：}} \\\\ \n".$part4_2."{\\sym{位点描述：}} \\\\ \n".$part4_3;
}

################### fb #################
my ($fb1,$fb2,$fb3)=("","","");
$part6_0='';
if(%site_dis_fb)
{
#	$fb_num++ if($negative != 2);
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
			($c1,$c2,$c3,$c4,$c5,$c6,$c6_2,$s2)=split(/\t/, $gene_site{$key1}{$key2});
			my $c11=word_wrap($c1,10);
			my $c22=word_wrap($c2,12);
			my $c33=word_wrap($c3,12);
						
			my ($c7,$c8,$c9)=('','','');
			my ($s3_1,$s3_2)=('','');
			my $part2_2='';
			my $part2_flag;
			my $j=0;
			
			foreach my $key3 (sort { $site_dis_fb{$key1}{$key2}->{$a} <=> $site_dis_fb{$key1}{$key2}->{$b} } keys %{ $site_dis_fb{$key1}{$key2} })
			{
				$part2_flag=0;
				my @arr=@{ $site_dis_fb{$key1}{$key2}{$key3} };
				$c8='';
				$c9=$key3;

				for(my $i=0;$i<@arr;$i++)
				{
					#next if($arr[$i][1] ne $key3);
					if($i>0 && $arr[$i][3] ne $arr[0][3]){ $part2_flag=1; }
					$c7   = $arr[$i][3];

					if($arr[$i][1] eq "-" && $arr[$i][0] ne "-")
					{
						$arr[$i][1]="不明";
					}
					$c8="$arr[$i][0]($arr[$i][1])" if($i==0);
					$c8.="\\\\$arr[$i][0]($arr[$i][1])"             if($i>0 && $part2_flag==0);
					$c8.="\\\\$arr[$i][0]($arr[$i][1],$arr[$i][3])" if($i>0 && $part2_flag==1);
					$c9=$key3; ###$arr[$i][4];			
	
					#$arr[$i][3]=~s/\~/--/g if($arr[$i][3]=~/\~/);
				}
			
				my $max_row=keys %{$site_dis_fb{$key1}{$key2}};
				#$s3.="$c1, $c2, $c3 \\\\ \n".$s3_1.": ".$s3_2;
			
				#$part2_2=qq[ & $c7 & \\makecell[{}{m{34mm}}]{$c8} & \\makecell[{}{m{9.5mm}}]{$c9} \\\\ \\cline{7-9}\n ] if($j==0 && $max_row>1);
			
				$part2_2=qq[ & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8}  \\\\ \\cline{7-9}\n ] if($j==0 && $max_row>1);
				$part2_2=qq[ & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8}  \\\\ \\hline\n ] if($j==0 && $max_row==1);
				$part2_2.=qq[ & & & & & & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8} \\\\ \\cline{7-9}\n ] if($j>0 && $j+1!=$max_row);
				$part2_2.=qq[ & & & & & & $c9 & $c7 & \\makecell[{}{m{34mm}}]{$c8} \\\\ \\hline\n ] if($j>0 && $j+1==$max_row);
			
				$j++;
			}
		
			if($j==1){
				$fb1.=qq[ $c11 & $c22 & $c33 & \\makecell[c]{$c4} & $c5 & \\makecell[{}{m{5.8mm}}]{$c6} ].$part2_2;
			}else{
				$fb1.=qq[ \\multirow{$j}{*}{$c11} & \\multirow{$j}{*}{$c22} & \\multirow{$j}{*}{$c33} & \\multirow{$j}{*}{\\makecell[c]{$c4} } & \\multirow{$j}{*}{$c5} & \\multirow{$j}{*}{ \\makecell[{}{m{5.8mm}}]{$c6} } ].$part2_2;
			}
		
		}
		
=head		
			@arr=split(/\t/, $gene_site{$key1}{$key2});
			$arr[0]=word_wrap($arr[0],10);
			$arr[1]=word_wrap($arr[1],12);
			$arr[2]=word_wrap($arr[2],12);
			@arr2=@{$site_dis_fb{$key1}{$key2}};
			my $i=0;
			my ($acmg,$zbx);
			for($i=0;$i<@arr2;$i++)
			{
				$arr2[$i][0]=~s/\(/\\\\\(/g;
				$arr2[$i][1]="不明" if($arr2[$i][1] eq "-" && $arr2[$i][0] ne "-");
				$fb2="$arr2[$i][0]($arr2[$i][1])" if($i==0);
				$fb2.="\\\\$arr2[$i][0]($arr2[$i][1])" if($i>0 and $fb2!~$arr2[$i][0]);#wph del 181010
				$acmg=$arr2[$i][2];
				$acmg=~s/,/, /g;
				$zbx =$arr2[$i][3];
			}
			$fb1.=qq[ $arr[0] & $arr[1] & $arr[2] & \\makecell*[c]{$arr[3]} & $arr[4] & $arr[5] & $acmg & $zbx & \\makecell[{}{m{34mm}}]{$fb2} \\\\ \\hline\n ];
=cut		
		
		$fb3.="{\\sye{$key1：}} \\\\ \n"."{\\sym{基因描述：}} \\\\ \n"."$gene_fb{$key1} \\\\ \n"."{\\sym{疾病描述：}} \\\\ \n";
		foreach my $key3 ( keys %{$disease_fb{$key1}} )
		{
			$fb3.=$disease_fb{$key1}{$key3}." \\\\ \\\\ \n";
		}
		$fb3.="\\\\ \n";	
	}
	$part6_0="\\newpage \n \\begin{spacing}{1.2} \n {\\sym{三、次要发现：}} \\\\ \n";
	$part6_0.=<<'END';
下表为无ACMG证据、证据不足、证据矛盾（偏致病类证据与偏良性证据共存）或偏良性证据的位点（可能与相关疾病出现的表型有关，具体情况请结合临床综合分析）\\
\end{spacing}
\fontsize{7}{8.4}\selectfont
\begin{longtable}{| C{13mm} | C{15mm} | C{15mm} | C{15mm} | C{6mm} | C{6mm} | m{10mm}<{\centering} | C{16mm} | m{34mm} | } \hline
{\sym{变异基因}} & {\sym{核苷酸变异}} & {\sym{氨基酸变异}} & \makecell[c]{ \sym{转录本}\\ \sym{外显子} } & {\sym{变异状态}} & {\sym{变异类型}} & \makecell[c]{ \sym{ACMG} \\ \sym{条目} } & {\sym{致病性}} & \makecell[c]{\sym{相关疾病(遗传模式)}} \\ \hline
END
#wph add makecell[{}{m{34mm}}] 181212
	$part6_0.=$fb1."\\end{longtable} \n";
	#print $part6_0;
}

$part6_0.=<<'END';
%\vspace*{3mm}
\begin{spacing}{1.3}
\zihao{6} \color{MyFontGray}
{\sym{备注：}} \\
1. {\sym{遗传模式：}}疾病对应的遗传模式，分为 AD（常染色体显性遗传）、AR（常染色体隐性遗传）、XL(X 染色体连锁)、XD/XLD（X 染色体显性遗传）、 XR/XLR（X 染色体隐性遗传）。\\
2. {\sym{核苷酸变异：}}\\
\hspace*{3mm} 1) 点变异的核酸改变表示： c.**（数字）碱基>碱基，或 c.**（数字）+**（数字）碱基>碱基。如：\\
\hspace*{7mm} a) c.387G>A，表示：该基因的编码序列第 387位的碱基 G变异为 A；\\
\hspace*{7mm} b) c.510+9T>-，表示：该基因的编码序列第 510位的碱基后的第 9 个碱基 T（位于内含子区域）缺失，距离该外显子边界（c.510）9bp 的位置。\\
\hspace*{3mm} 2) 插入缺失变异（Indel）位置的核酸改变表示：c.**(起始位置数字)\_(终止位置数字)：插入/缺失&碱基。如：c.4030\_4032delGAG，表示：从编码序列第 4030 位碱基到 4032位碱基位置，缺失3个碱基GAG 。\\
3. {\sym{氨基酸改变：}}核苷酸变异会导致其所编码的氨基酸产生变化，具体分为：\\
\hspace*{3mm}   1) {\sym{错义变异：}}碱基的变异导致其所在氨基酸的改变。如：p.氨基酸（数字）氨基酸，如 p.Met129Ile，表示：该基因所编码的第 129 个氨基酸由Met(蛋氨酸)变为 Ile(异亮氨酸)；\\
\hspace*{3mm}   2) {\sym{无义变异：}}指由于某碱基的变化导致其所在氨基酸变成终止密码子。表示为：p.氨基酸(数字) * 。如p.Gln116*，表示：该基因所编码的第116个氨基酸由 Gln（ 谷氨酰胺） 变为终止密码子；\\
\hspace*{3mm}   3) {\sym{剪接位点变异：}}这种变异会影响基因的转录过程，导致转录产物 mRNA 序列的异常，从而影响蛋白功能。这类变异不会直接导致氨基酸的变化，故而标记为 ` $-$ ' ， 一般位于剪切位点±1/2的变异，对剪切影响较大。\\
\hspace*{3mm}   4) {\sym{移码变异：}}指由于插入或缺失 N 个（非 3 整数倍）碱基的变异，这类变异会导致读码框发生错位，从而导致该位置及之后的氨基酸序列的明显改变。表示为：p.氨基酸（数字1）氨基酸 fs * (数字2)，fs 表示 frameshift（移码）。数字 2 表示从突变位置开始算起，后面还可以翻译但框移的氨基酸的数量。如p.Arg311LysfsTer7，即表示：从第 311 号氨基酸 Arg 开始的氨基酸合成发生改变，并在改变后的第 7 个氨基酸终止; \\
\end{spacing}
END

#$part6_0.= "\\\\ \\\\ \\\\ \n \\begin{spacing}{1.3}  \\fontsize{10.5}{12.6}\\selectfont \n".$fb3."\\end{spacing}";
$part6_0.= "\\\\ \n";

################### IDT_stat #################
if($types eq "WES")
{
	open IN, "$genelist" or die $!;
	my $i=0;
	my $i_num=0;
	my $step=0;
	my @cut_gene;
	my %tmp_gene;
	while(<IN>)
	{
		chomp;
		@cut=split /\t/,$_;
		if($cut[0] eq $samples )
		{
		#	$disease_name=$cut[1] if($i==0 && $cut[1] ne "高血压");
		#	$disease_name.=",$cut[1]" if($i>0 && $cut[1] ne "高血压");
			@cut_gene=split /, /,$cut[3];
			foreach my $key (@cut_gene)
			{
				next if($key eq "-");
				if($i_num%100==0 && int($i_num/100)>$step && $i_num>0){
					$step=$i_num/100;
					$disease_gene_line=~s/; $//;
					$disease_gene_line_hash{$step}=$disease_gene_line."等";
					$disease_gene_line="";
				}
				if(!exists $tmp_gene{$key})
				{
					$disease_gene_line.="$key; ";
					$tmp_gene{$key}=1;
					$i_num++;
				}

			}
			$i++;
		}
		#$disease_name=~s/单基因糖尿病,//g;
	}
	if($disease_gene_line ne "")
	{
		$step++;
		$disease_gene_line=~s/; $//;
		$disease_gene_line_hash{$step}=$disease_gene_line."等" if($step >1);
		$disease_gene_line_hash{$step}=$disease_gene_line      if($step==1);
	}
	$disease_name=$i;
	$disease_gene=keys %tmp_gene;
	close IN;
}

################### negative #################
if($negative==2)
{
	$part1=<<'END';
\noindent \\ \\ %\\
\begin{spacing}{1.3}
{\sym{二、主要变异检测结果：}} \\
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
		$part1.="此次全外显子基因检测，包含".$disease_name."个疾病相关的".$disease_gene."个基因。";
	}
	$part1.=<<'END';
检测结果显示您并未发生有风险的位点变异，因此，您的疾病症状并不是由这些基因的变异引起的，可能是目前研究尚未涉及的基因变异，也可能是由于其他因素，如不良生活习惯、环境、其他疾病等引起。建议您注意生活细节，积极改变不良生活方式，并配合医生治疗。\\ \begin{spacing}{0.8} \end{spacing}
以上检测结果仅基于当前的科学研究，随研究文献的更新，疾病-基因的对应关系可能会发生变化，因此我们的报告具有一定的时效性。\\
\end{spacing}
\noindent
END

}
	
if($negative==1)
{
	$part1=<<'END';
\noindent \\ \\ %\\
\begin{spacing}{1.3}
{\sym{二、主要变异检测结果：}} \\
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
检测结果显示您并未发生高风险的位点变异；一些变异位点由于相关研究较少、ACMG证据条目不足，基于现有研究判断其致病性等级较低，但不排除导致疾病的可能性，这类位点的致病性等级及相关证据条目展示在附表中，供您和临床医生参考。\\ \begin{spacing}{0.8} \end{spacing}
若检出位点的相关疾病与您现有的临床表型不符，那么您的相关临床症状可能与其他未知基因变异、疾病的不完全外显性、不良生活习惯或环境等因素相关。建议您注意生活细节，积极改变不良生活方式，或去医院进行相关疾病的进一步检测，听从医生的临床建议，定期体检，并采取积极的生活干预。\\ \begin{spacing}{0.8} \end{spacing}
若您未发现任何相关疾病表型，则这类 ACMG 致病类证据较少的位点可能并不致病。\\ \begin{spacing}{0.8} \end{spacing}
以上检测结果仅基于当前的科学研究，随研究文献的更新，位点致病性情况可能会发生变化，属于正常现象。如果您想进一步了解基因变异的信息或家族遗传史，建议联系家人一起进行检测。\\ \begin{spacing}{0.8} \end{spacing}
\end{spacing}
\noindent
END
	
}

################### list #################
my %sample_uniq;
if($negative==1 or $negative==0 or $negative==2 or $negative==3)
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
\makecell[{{p{8.1cm}}}]{新生儿糖尿病(NDM)\\Neonatal Diabetes Mellitus} & {KCNJ11 、INS 、GCK 、PDX1 、MNX1 、GLIS3 、NEUROD1 、NEUROG3 、PAX6 、PTF1A 、RFX6 、GATA4 、NKX2-2 、SLC19A2 、IER3IP1 、ZFP57 、ABCC8 、HNF1B} \\ \hline
\makecell[{{p{8.1cm}}}]{其他单基因糖尿病} & {SLC2A2 、WFS1 、TRMT10A、PPP1R15B 、EIF2AK3 、CISD2 、INSR 、ALMS1 、BLM 、HADH 、FOXP3、FXN} \\ \hline
\end{longtable}
\noindent \\
END
	}#IER31P1 、PLAGL1 wph 20190423删除 TRMT10A、PPP1R15B拆分  sxq  add FXN 20190425
	
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
			#$part1_list.="\\makecell[{{p{8.2cm}}}]{$cut[0]} & \\makecell[{{p{7cm}}}]{$cut[1]} \\\\ \\hline \n";
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
		$part1_list.="{\\sym\\sye{附表 1 ：}}";
		$part1_list.=<<'END';
本次检测疾病相关基因如下\\
\renewcommand\arraystretch{1.5}
\fontsize{9}{10.8}\selectfont
\begin{longtable}{ L{0.1cm}  L{15.2cm} }
\hline
%\makecell[c]{\sym{序号}} & \makecell[c]{\sym{相关基因}} \\ \hline
\makecell[c]{\sym{}} & \makecell[c]{\sym{相关基因}} \\ \hline
END

		open IN, "$genelist" or die $!;
		while(<IN>)
		{
			chomp;
			@cut=split /\t/,$_;
			
			if($cut[0] eq $samples )
			{
				if ($cut[1]=~/心源性猝死/){$part1_list.="\\makecell[{{p{4cm}}}]{心源性猝死} & \\makecell[{{p{11cm}}}]{ABCA3, ABCC9, ACTA1, ACTA2, ACTC1, ACTN2, ACVRL1, AKAP9, ALDH1A2, ALG10, ANK2, ANKRD1, APOB, BAG3, BMPR2, BRAF, BVES, CACNA1C, CACNA2D1, CACNA2D4, CACNB2, CALM1, CALM2, CALM3, EIF2AK4...} \\\\ \\hline \n";next;}
				if ($cut[1]=~/单基因高血压|单基因糖尿病/){next;}
				if ($cut[1] eq "高血压") {
					$cut[1] = "单基因高血压";
					$part1_list.="\\makecell[{{p{4cm}}}]{单基因高血压} & \\makecell[{{p{11cm}}}]{AIP, ARMC5, BSND, CACNA1D, CACNA1H, CASR, CDKN1B, CLCN2, CLCNKA, CLCNKB, CUL3, CYP11B1, CYP17A1, EGLN1, EGLN2...} \\\\ \\hline \n";
					next;
				}
				if(!exists $sample_uniq{$samples} || $sample_uniq{$samples}!=1)
				{
					foreach my $key ( sort {$a<=>$b}  keys %disease_gene_line_hash)
					{
						#$part1_list.="\\makecell[{}{p{4cm}}]{$cut[1]} & \\makecell[{{p{11cm}}}]{$cut[2]} \\\\ \\hline \n";
						#$part1_list.="\\makecell[{}{p{1cm}}]{$key} & \\makecell[{{p{14cm}}}]{$disease_gene_line_hash{$key}} \\\\ \\hline \n";
						$part1_list.="\\makecell[{}{p{0.1cm}}]{} & \\makecell[{{p{15cm}}}]{$disease_gene_line_hash{$key}} \\\\ \\hline \n" if($key == 1);
					}
					$sample_uniq{$samples}=1;
				}
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
$part6="\\fontsize{10.5}{12.6}\\selectfont\n{\\sym\\sye{附表 2 ：}} ";
$part6.=<<'END';
本次检测的主要测序参数如下\\
本次检测采用方法为目标序列捕获高通量测序技术，采用 Illumina 公司 NovaSeq 系列测序平台，产出数据的质量及准确性皆为业内最高标准，配合严谨的信息分析流程和数据解读方法，确保您宝贵的基因信息得到最科学的管理。\\
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

$check="";

#######################################################
sub word_wrap{ #wph add 190628
	my $word=shift;
	my $linewid=shift;
	my $lenword=length($word);
	my $flag=0;
	if($word=~/;/){
		$word=~s/;/;\\\\/g;
	}elsif($lenword>$linewid && $lenword<2*$linewid){
		if($word =~/\_/ && $word =~/del/){
			if(index($word,"del")>$linewid){
				$flag=1;
			}else{
				$word=~s/del/\\\\del/g;
			}
		}elsif($word =~/\_/ && $word =~/ins/){
			if(index($word,"ins")>$linewid){
				$flag=1;
			}else{
				$word=~s/ins/\\\\ins/g;
			}
		}elsif($word =~/delins/ ){
			$word=~s/delins/delins\\\\/g;
		}elsif($word =~/fs/ ){
			$word=~s/fs/\\\\fs/g;
		}elsif($lenword>$linewid){
			$word=substr($word,0,$linewid)."\\\\".substr($word,$linewid,$lenword-$linewid);
		}
	}
	if($lenword>=2*$linewid or $flag){
		my $cycleall=$lenword/$linewid;
		my $siteindex=$linewid;
		if($word=~/ins/){
			$siteindex=index($word,"ins");
		}elsif($word=~/del/){
			$siteindex=index($word,"del");
		}
		my @replword;
		foreach my $cycle((0..$cycleall)){
			if($siteindex<$linewid){
				if($cycle==0){
					push @replword,substr($word,0,$siteindex);
				}else{
					push @replword,substr($word,($cycle-1)*($linewid-1)+$siteindex,$linewid-1);
				}
			}else{
				push @replword,substr($word,$cycle*($linewid-1),$linewid-1);
			}
		}
		push @replword,substr($word,$cycleall*$linewid,$lenword-$cycleall*$linewid);
		$word=join("\\\\",@replword);
	}
	if($word=~/\\\\/){
		$word="\\makecell[c]{$word}";
	}
	
	return $word;
}


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
                if($date=~/2019心康送样信息单-(\d+)\.xlsx/)
                {
                        push @all_date,$1;
                }
        }

        $lastest_date= max(@all_date);
        my $lastest_sitedb_path ="$product_site/2019心康送样信息单-$lastest_date.xlsx";
	my $lastest_sitedb_path2="$product_site/2019心康送样信息单-$lastest_date.txt";
        print "$lastest_sitedb_path\n";
        return $lastest_sitedb_path,$lastest_sitedb_path2;
}

sub judge_lastest_sitedb2{
	my $product_type=shift;
	my %hash_product_sitepath=();
	my @all_date;
	my $lastest_date;
	my $lastest_sitedb_path;

    my $gene_descript_path="/PUBLIC/pipline/database/Auto_pipline/gene_description/";#wangpenghui
    my $disease_descript_path="/PUBLIC/pipline/database/Auto_pipline/disease_description/";#wangpenghui

    $hash_product_sitepath{'Gene'}=$gene_descript_path;
    $hash_product_sitepath{'Disease'}=$disease_descript_path;
	if( exists $hash_product_sitepath{$product_type}){
		my @all_path=glob "$hash_product_sitepath{$product_type}/*";
		foreach(@all_path){
			my $date=(split/\//)[-1];
			if($date=~/^\d+$/i){
				push @all_date,$date;
			}
		}
		$lastest_date=max(@all_date);
		my @lastest_sitedb_path=glob "$hash_product_sitepath{$product_type}/$lastest_date/*right.txt";
		foreach(@lastest_sitedb_path){
			if(/right/i){
				$lastest_sitedb_path=$_;
			}
		}
	}
	print "$lastest_sitedb_path\n";
	print "$product_type\n";
	return $lastest_sitedb_path;
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
        #$biaoti="单基因遗传病基因检测报告";
        $biaoti="全外显子基因检测报告";
}
$biaoti2=qq[\\section*{\\sym{$biaoti}}
\\zihao{-4} {\\sym{检测编号：$samples}} \\\\ \\\\
\\color{MyFontGray} \\fontsize{10.5}{12.6}\\selectfont
{\\sym{一、样本信息：}} \\\\ \\\\
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
	$part7_1="1. 检测结果仅对送检的该样本负责。\\\\ \n2. 本次检测重点关注疾病数据库收录的与".$disease_name."疾病相关的".$disease_gene."个基因的全部外显子区域";
}


$part7=<<'END';
\newpage
\vspace*{2mm}
\begin{spacing}{1.4}
\zihao{-4}{\sym {\underline {进一步检测建议：} } }\\ \\
\fontsize{9}{10.8}\selectfont
如果送检医生认为本报告内容不足以解释临床表型，可考虑进行进一步检测：\\
1. 线粒体测序：全外显子检测无法检测到线粒体基因变异。如考虑到线粒体基因变异致病的可能性，可进一步进行线粒体基因检测。\\
2. CNV 检测：全外显子检测无法非常准确、特异地检测大片段变异（ CNV）。如考虑大片段变异的可能性，可进一步进行 CNV 相关检测。\\
3. 表观遗传学检测：外显子变异无法展示由于甲基化、乙酰化等表观遗传学改变对基因及基因产物的影响，可进一步进行表观遗传学检测。\\
\end{spacing}
\vspace*{3mm}
\begin{spacing}{1.4} 
\zihao{-4}{\sym {\underline {温馨提示：} } }\\ \\
\fontsize{9}{10.8}\selectfont
END
$part7.=$part7_1."、外显子与内含子衔接处的片段，不包括基因调控区及深度内含子区，检测仅能提供上述覆盖区域内的位点变异和小片段插入/缺失（<25bp）等信息分析。此外，由于二代测序的局限性，现有的技术难以捕获部分GC含量高的基因区域，GC含量高的区域也会有未能覆盖的情况。\\\\ \n";
$part7.=<<'END';
3. 本检测不适用于分析基因组结构变异（如大片段拷贝数缺失、插入或倒位）、动态变异、复杂重组等特殊类型的变异。MLPA 技术用于检测基因的大片段变异。\\
4. 以上结论均为实验室检测数据，用于疾病筛查之目的，仅供参考，不作为最终诊断结果，相关解释请咨询临床医生。\\
5. 本检验所对该结果保密并依法保护用户的隐私，如因受检者个人原因出现的信息外泄，本检验所不承担相应责任。\\
6. 本检验所保留对上述结果的最终解释权，如有疑义，请在收到报告后的 7 个工作日内与我们联系。\\
\end{spacing}
\vspace*{5mm}
END

my $all_tex;
if($negative == 0)
{
	$all_tex=$yemei.$biaoti2.$info_sample.$part1.$part2."\n \\end{longtable}\n".$part3."\\fontsize{10.5}{12.6}\\selectfont \n%\\noindent \\\\ \n\\begin{spacing}{1.3}\n{\\sym {检测结果说明：}} \\\\ \n".$part4."\\end{spacing}\n".$part5.$part6_0.$part9.$part10.$part7.$part1_list.$part6.$part8.$check."\n\\end{document} \n";
}
elsif($negative == 1)
{
	$all_tex=$yemei.$biaoti2.$info_sample.$part1.$part6_0.$part9.$part10.$part7.$part1_list.$part6.$part8.$check."\n\\end{document} \n";
}
elsif($negative == 2)
{
	$all_tex=$yemei.$biaoti2.$info_sample.$part1.$part2.$part9.$part10.$part7.$part1_list.$part6.$part8.$check."\n\\end{document} \n";
}
elsif($negative == 3)
{
	$all_tex=$yemei.$biaoti2.$info_sample.$part1.$part2."\n \\end{longtable}\n".$part3."\\fontsize{10.5}{12.6}\\selectfont \n%\\noindent \\\\ \n\\begin{spacing}{1.3}\n{\\sym {检测结果说明：}} \\\\ \n".$part4."\\end{spacing}\n".$part5.$part9.$part10.$part7.$part1_list.$part6.$part8.$check."\n\\end{document} \n";
}

open OUT, ">$out/$samples.$types.tex";
print OUT $all_tex;
close OUT;
