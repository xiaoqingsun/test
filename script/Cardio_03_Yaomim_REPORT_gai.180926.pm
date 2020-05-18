#!/usr/bin/env perl
#===============================================================================
#
#         FILE: Report.GeneID.pl
#
#        USAGE: ./Report.GeneID.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (),
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 2016年10月25日 12时17分51秒
#     REVISION: ---
#===============================================================================
package Cardio_03_Yaomim_REPORT;
use strict;
#use warnings;
use utf8;
use Encode;
use POSIX;
use Data::Dumper;

my (%hash,%hash1,%hash2,%hash3,%hash4,%hash_only,%hashyn);
#my (%hash,%hash2,%hash3,%hash4,%hash_only,%hashh,%hasht,$references);
sub report{
	my ($gresult,$yresult,$dresult,$part)=@_;
############################### main ##############################
	my $tex=&shuru($gresult,$yresult,$dresult);
	return $tex;
#################读入文件进行处理##############

sub shuru{
	my ($file,$file2,$file3)=@_;
	my @path=split/\//,$file;
    my @sample_path=split/\./,$path[-1];
    my $sample=$sample_path[0];
	open IN1,"<:encoding(utf-8)","$file" || die "can't open the $file!";      #高血压用药部分PDF文件
	while(<IN1>){
		chomp;
		my $level;
		my @cut=split/\t/,$_;
		next if(/^药物类别/ || /^#/);
		if($cut[0]=~/$sample/){
			my $tttt=shift @cut;                                                  #针对第一列多了#样本编号
		}
		next if($cut[4]=~/\-/ or $cut[4]=~/正常用药/);
		next if($cut[6] eq "-");
		my @ref=split/\|/,$cut[11];
		$hashyn{$cut[2]}{$ref[0]}++;
		$hashyn{$cut[2]}{$ref[1]}++;
		if($hashyn{$cut[2]}{$ref[0]} <2){
			push @{$hash1{$cut[2]}},$ref[0];		             #处理参考文献的
		}
		if($hashyn{$cut[2]}{$ref[1]} <2){
			push @{$hash1{$cut[2]}},$ref[1];
		}
		my @a=split//,$cut[9];
        my $gety=join "/",@a;
		if($cut[4]=~/不建议使用/){
			$level="不建议使用";
		}
		elsif($cut[4]=~/增加用药、小心使用/){
			$level="小心使用";
		}
		elsif($cut[4]=~/增加用药/){
			$level="增加剂量";
		}
		elsif($cut[4]=~/小心使用/){
			$level="小心使用";
		}
		my @val=($cut[2],$level,$cut[3],$cut[5]);     #药物、建议等级、其他名称、用药建议
		my $key=join "\t",@val;
		$hash_only{$key}++;
		if($hash_only{$key} < 2){
			push @{$hash3{$cut[0]}},[@val];
		}
		push @{$hash{$cut[0]}{$key}},$cut[6];        #一维数组，值为结果分析
		my @vall=($cut[7],$cut[8],$gety,$cut[10]);        #基因、位点（rs_number）、基因型、位点结论
		push @{$hash2{$cut[2]}},[@vall];     #已药物名称做键，获得它相关的基因、位点信息
	}
	open IN2,"<:encoding(utf-8)","$file2" || die "can't open the $file2!";    #华法林与氯吡格雷药物文件
    while(<IN2>){
		chomp;
		my $level;
		my @cut=split/\t/,$_;
		next if/^sample/;
		next if($cut[12]=~/正常使用/);
		$cut[8]=~s/心血管类/心血管病常用药/g;
		my @ref=split/\|/,$cut[16];
		$hashyn{$cut[9]}{$ref[0]}++;
		$hashyn{$cut[9]}{$ref[1]}++;
        if($hashyn{$cut[9]}{$ref[0]} <2 or $hashyn{$cut[9]}{$ref[1]} <2){
            push @{$hash1{$cut[9]}},$ref[0];                    #处理参考文献的
			push @{$hash1{$cut[9]}},$ref[1];
        }
		my @a=split//,$cut[11];
		my $gety=join "/",@a;
		if($cut[12] eq "建议换药"){
			$level="不建议使用";
		}
		elsif($cut[12] eq "调整剂量"){
			$level="调整剂量";
		}
		elsif($cut[12] eq "酌情使用"){
			$level="小心使用";
		}
		my @val=($cut[9],$level,$cut[10],$cut[14]);
		my $key=join "\t",@val;
		$hash_only{$key}++;
		if($hash_only{$key} < 2){
			push @{$hash3{$cut[8]}},[@val];
		}
		push @{$hash{$cut[8]}{$key}},$cut[15];
		my @vall=($cut[6],$cut[7],$gety,$cut[19]);
		push @{$hash2{$cut[9]}},[@vall];
	}
	open IN3,"<:encoding(utf-8)","$file3" || die "can't open the $file2!";    #新增药物 2017.11.24
    while(<IN3>){
		chomp;
		my $level;
		my @cut=split/\t/,$_;
		next if/^sample/;
		next if($cut[12]=~/正常使用/);
		$cut[8]=~s/心血管类/心血管病常用药/g;
		my @ref=split/\|/,$cut[16];
		$hashyn{$cut[9]}{$ref[0]}++;
		$hashyn{$cut[9]}{$ref[1]}++;
        if($hashyn{$cut[9]}{$ref[0]} <2 or $hashyn{$cut[9]}{$ref[1]} <2){
            push @{$hash1{$cut[9]}},$ref[0];                    #处理参考文献的
			push @{$hash1{$cut[9]}},$ref[1];
        }
		my @a=split//,$cut[11];
		my $gety=join "/",@a;
		$cut[12]=~s/\s//g;
		if($cut[12] eq "建议换药"){
			$level="不建议使用";
		}
		elsif($cut[12] eq "调整剂量"){
			$level="调整剂量";
		}
		elsif($cut[12] eq "酌情使用"){
			$level="小心使用";
		}
		my @val=($cut[9],$level,$cut[10],$cut[14]);
		my $key=join "\t",@val;
		$hash_only{$key}++;
		if($hash_only{$key} < 2){
			push @{$hash3{$cut[8]}},[@val];
		}
		push @{$hash{$cut[8]}{$key}},$cut[15];
		my @vall=($cut[6],$cut[7],$gety,$cut[20]);
		push @{$hash2{$cut[9]}},[@vall];
	}
	my @med_class=("利尿剂","β受体阻断剂","血管紧张素转换酶抑制剂ACEI","血管紧张素II受体拮抗剂ARB","钙通道阻断剂","高血糖","他汀类药物","心血管类","高血脂","心血管病常用药");
	
	#my @med_class=("利尿剂","β受体阻断剂","血管紧张素转换酶抑制剂ACEI","血管紧张素II受体拮抗剂ARB","钙通道阻断剂","心血管病常用抗凝药","双胍类药物","磷酸脲类促胰岛素分泌药","胰岛素增敏类药物","肠促胰岛素类似物","他汀类药物","贝特类药物","胆固醇吸收抑制剂");
	close IN1;
	close IN2;
	close IN3;
my $tail=qq[
\\end{YMpart}
];
my $tailr=<<'END';
\pagestyle{fancy}
%\newsavebox{\mygraphic}
%\arrayrulecolor{white}
%\sbox{\mygraphic}{\mbox{
%\begin{tabular}{L{12.5cm} L{4cm}}
%\cellcolor{LightBlue}\color{DarkBlue}{\fontsize{23}{27.6}{\sym{心血管用药指导科普}}} &
%\multirow{2}*{\includegraphics[height=3.75cm,width=4cm]{7-3.pdf}} \\ [1.8cm]\cmidrule{1-1}
%\cellcolor{LightBlue}\color{MyFontGray}{\fontsize{15}{18}{高血压药物基因检测}} & \\ [0.77cm]
%\end{tabular}}}
%\fancyhead[C]{\usebox{\mygraphic}}
%\fancyfoot[RO,LE]{\zihao{6} Novo\_\thepage}
\noindent
%\fontsize{6.5}{7.8}{\color{MyFontGray}{
\setlength{\columnsep}{2pc}
\begin{multicols}{2}
\zihao{-6}
END
	my $report=&REPORT(@med_class);
	my $Latex="$report"."$tail";
	return $Latex;
}
##################子程序################

sub REPORT{
	my @med_class=@_;
	my $tabinfo;

my $title=<<'END';
%\arrayrulecolor{DarkBlue} \hline
\vspace*{0.5cm}
\begin{spacing}{1.5}
\fontsize{12}{14.4}\selectfont
\color{MyFontGray}
\noindent
END

my $head=qq[
\\newpage
\\addcontentsline{toc}{subsection}{ⅱ 检测结果详细解析}
\\begin{center}
\\zihao{3}{\\color{DarkBlue} \\sye \\sym 检测结果详细解析}
\\end{center}
];
	foreach my $a(@med_class){
		if (exists $hash3{$a}){
			my @co=@{$hash3{$a}};
			foreach my $b(@co){
				my @arr=@$b;
				my $references;
				foreach(@{$hash1{$arr[0]}}){            #处理参考文献，药物作为键
					my $cankao=qq[
$_ \\par ];
					$references.="$cankao";
				}
				my $k=join "\t",@arr;
				my @value=@{$hash{$a}{$k}};     #获得它的值即：结果分析的数组
				my @c = grep{++$hash4{$_}<2} @value;
				my $analyse=join ".",@c;
				my $result=&TABLE($arr[0]);
				my $color;
				my $ym;
				if($arr[1] =~ /不建议使用/){
			        $color='red';
					$ym='建议换药';
			    }
			    elsif($arr[1] =~ /增加剂量/){
				    $color='orange';
					$ym='调整剂量';
			    }
				elsif($arr[1] =~ /调整剂量/){
					$color='orange';
					$ym='调整剂量';
				}
			    elsif($arr[1] =~ /小心使用/){
			        $color='yellow';
					$ym='小心使用';
			    }
				my $tabin=qq[\\textcolor{DarkBlue}{\\sym{药物：}}$arr[0]\\\\
\\textcolor{DarkBlue}{\\sym{建议等级：}}\\includegraphics[height=0.3cm]{药敏$ym.pdf}\\textcolor{$color}{$arr[1]}\\\\
\\textcolor{DarkBlue}{\\sym{其他名称：}}$arr[2]\\\\
\\textcolor{DarkBlue}{\\sym{用药建议：}}$arr[3]\\\\
\\end{spacing}
\\arrayrulecolor{DarkBlue} \\hline
\\vspace*{1.1cm}
];
my $tabin1=qq[
\\vspace*{-8mm}
\\noindent
\\color{MyFontGray}
\\fontsize{12}{14.4}{\\textcolor{DarkBlue}{\\sym{结果分析：}}}\\\\ \\\\
\\fontsize{9.5}{11.4}\\selectfont
$analyse  \\par];
my $reftex;
if($references eq ""){
	$reftex="";}
else{
	$reftex=qq[
\\vspace*{6mm}
\\fontsize{12}{14.4}{\\textcolor{DarkBlue}{\\sym{参考文献：}}} \\par
\\fontsize{9}{10.8}\\selectfont
$references];}
				$tabinfo.="$title"."$tabin"."$result"."$tabin1"."$reftex"."\\newpage\n";
			}
		}
	}
	my $tableinfo="$head"."$tabinfo";
	return $tableinfo;
}

sub TABLE {
	my(%hashh,%hasht);
	my $m=shift;
	my @med=@{$hash2{$m}};    #将药物对应的所有基因/位点等放到数组@med
	my $i;
my $td=<<'END';
\\
\noindent
\fontsize{12}{14.4}{\textcolor{DarkBlue}{\sym{检测结果：}}}\\
\vspace*{-0.5cm}
\begin{center}
\fontsize{9}{10.8}\selectfont
\color{MyFontGray}
\arrayrulecolor{DarkBlue}
\renewcommand\arraystretch{1.5}
\begin{VerbatimOut}{\jobname.vrb}
\begin{longtable}{ p{.15\textwidth} p{.15\textwidth} p{.2\textwidth}<{\centering} p{.4\textwidth} }
\rowcolor{DarkBlue}
\fontsize{9.5}{10.4}{\textcolor{white}{\sym{基因}}} & \fontsize{9.5}{10.4}{\textcolor{white}{\sym{位点}}} & \fontsize{9.5}{10.4}{\textcolor{white}{\sym{您的基因型}}} & \fontsize{9.5}{10.4}{\textcolor{white}{\sym{位点结论}}}  \tabularnewline
END

my $tail=<<'END';
\end{longtable}
\end{VerbatimOut}
\LTXtable{\textwidth}{\jobname.vrb}
\end{center}
END
	my ($tabin2,$tableinfo,$info);
	for($i= 0;$i<=$#med;$i++){
		$hashh{$med[$i]->[0]}++;
		$hasht{$med[$i]->[3]}++;
	}
	my $h=$hashh{$med[0]->[0]};                      #针对一个药多位点多基因情况
	my $t=$hasht{$med[0]->[3]};
	if($h==@med){
		if($t==@med){
			for($i= 0;$i<=$#med;$i++){
	  			if($i==0 and $t!=1){
					$tabin2=qq[\\multirow{$h}{*}{\\cellcolor{white} $med[$i]->[0]} & \\cellcolor{white} $med[$i]->[1] & \\cellcolor{white} $med[$i]->[2] & \\multirow{$h}{*}{\\cellcolor{white} $med[$i]->[3]} \\tabularnewline \n];
					$tableinfo.="$tabin2";
				}
				if($i==0 and $t==1){
					$tabin2=qq[$med[$i]->[0] & $med[$i]->[1] & $med[$i]->[2] & $med[$i]->[3] \\tabularnewline \\hline \n];
                    $tableinfo.="$tabin2";
				}
				elsif(0<$i and $i<$#med and $i % 2 ==0){
					$tabin2=qq[ & \\cellcolor{white} $med[$i]->[1] & \\cellcolor{white} $med[$i]->[2] & \\tabularnewline \n];
					$tableinfo.="$tabin2";
				}
				elsif(0<$i and $i<$#med and $i % 2 ==1){
					$tabin2=qq[ & \\cellcolor{LightBlue} $med[$i]->[1] & \\cellcolor{LightBlue} $med[$i]->[2] & \\tabularnewline \n];
					$tableinfo.="$tabin2";
				}
				elsif($i==$#med and $i % 2 ==0){
					$tabin2=qq[ & \\cellcolor{white} $med[$i]->[1] & \\cellcolor{white} $med[$i]->[2] & \\tabularnewline \\hline \n];
					$tableinfo.="$tabin2";
				}
				elsif($i==$#med and $i % 2 ==1){
					$tabin2=qq[ & \\cellcolor{LightBlue} $med[$i]->[1] & \\cellcolor{LightBlue} $med[$i]->[2] & \\tabularnewline \\hline \n];
					$tableinfo.="$tabin2";
				}
			}
			$info="$td"."$tableinfo"."$tail";
		}
		else{
			for($i= 0;$i<=$#med;$i++){
				if($i==0){
					$tabin2=qq[\\multirow{$h}{*}{\\cellcolor{white} $med[$i]->[0]} & \\cellcolor{white} $med[$i]->[1] & \\cellcolor{white} $med[$i]->[2] & \\cellcolor{white} $med[$i]->[3] \\tabularnewline \n];
					$tableinfo.="$tabin2";
				}
				elsif(0<$i and $i<$#med and $i % 2 ==0){
					$tabin2=qq[ & \\cellcolor{white} $med[$i]->[1] & \\cellcolor{white} $med[$i]->[2] & \\cellcolor{white} $med[$i]->[3] \\tabularnewline \n];
					$tableinfo.="$tabin2";
				}
                elsif(0<$i and $i<$#med and $i % 2 ==1){
	                $tabin2=qq[ & \\cellcolor{LightBlue} $med[$i]->[1] & \\cellcolor{LightBlue} $med[$i]->[2] & \\cellcolor{LightBlue} $med[$i]->[3] \\tabularnewline \n];
					$tableinfo.="$tabin2";
				}
				elsif($i==$#med and $i % 2 ==0){
					$tabin2=qq[ & \\cellcolor{white} $med[$i]->[1] & \\cellcolor{white} $med[$i]->[2] & \\cellcolor{white} $med[$i]->[3] \\tabularnewline \\hline \n];
					$tableinfo.="$tabin2";
				}
				elsif($i==$#med and $i % 2 ==1){
					$tabin2=qq[ & \\cellcolor{LightBlue} $med[$i]->[1] & \\cellcolor{LightBlue} $med[$i]->[2] & \\cellcolor{LightBlue} $med[$i]->[3] \\tabularnewline \\hline \n];
					$tableinfo.="$tabin2";
				}
			}
		}
		$info="$td"."$tableinfo"."$tail";
	}
	else{
		for($i= 0;$i<=$#med;$i++){
			if($i<$#med and $i % 2 ==0){
				$tabin2=qq[$med[$i]->[0] & $med[$i]->[1] & $med[$i]->[2] & $med[$i]->[3] \\tabularnewline \n];
				$tableinfo.="$tabin2";
			}
			elsif($i<$#med and $i % 2 ==1){
				$tabin2=qq[\\rowcolor{LightBlue} $med[$i]->[0] & $med[$i]->[1] & $med[$i]->[2] & $med[$i]->[3] \\tabularnewline \n];
				$tableinfo.="$tabin2";
			}
			elsif($i==$#med and $i % 2 ==0){
				$tabin2=qq[$med[$i]->[0] & $med[$i]->[1] & $med[$i]->[2] & $med[$i]->[3] \\tabularnewline \\hline \n];
				$tableinfo.="$tabin2";
			}
			elsif($i==$#med and $i % 2 ==1){
				$tabin2=qq[\\rowcolor{LightBlue} $med[$i]->[0] & $med[$i]->[1] & $med[$i]->[2] & $med[$i]->[3] \\tabularnewline \\hline \n];
				$tableinfo.="$tabin2";
			}
		}
		$info="$td"."$tableinfo"."$tail";
	}
	return $info;
}

}
1;
