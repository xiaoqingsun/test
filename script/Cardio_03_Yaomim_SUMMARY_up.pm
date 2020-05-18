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
package Cardio_03_Yaomim_SUMMARY;
use strict;
use warnings;
use utf8;
use Encode;
use POSIX;
use Data::Dumper;


my %hash;

sub summary{
	my ($gresult,$yresult,$part,$module)=@_;
############################### main ##############################
    my $tex=&shuru($gresult,$yresult,$part,$module);
	return $tex;

#################读入文件进行处理##############
sub shuru{
	my $file=shift;
	my $file1=shift;
	my $part1=shift;
    my $module1=shift;
	my @path=split/\//,$file;
	my @sample_path=split/\./,$path[-1];
	my $sample=$sample_path[0];
	my $num=0;my $num1=0;my $num2=0;my $num3=0;
	my (%hash1,%hash2,%hash3);
	open IN1,"<:encoding(utf-8)","$file" || die "can't open the $file!";      #高血压用药部分PDF文件
	while(<IN1>){
		chomp;
		my $level;
		my @cut=split/\t/,$_;
        next if(/^药物类别/ || /^#/);
		if($cut[0]=~/$sample/){
			my $tttt=shift @cut;     #后来加的，针对第一列增加了#样本编号
		}
		next if($cut[4]=~/\-/);
		$hash2{$cut[2]}++;
		$cut[4]=~s/\s//g;
		if($hash2{$cut[2]}<=1 and $cut[4] eq "正常用药"){
			$level="正常使用";
		}
		elsif($hash2{$cut[2]}<=1 and $cut[4] eq "不建议使用"){
			$num1++;
			$level="不建议使用";
		}
		elsif($hash2{$cut[2]}<=1 and $cut[4] eq "增加用药、小心使用"){
			$num3++;
			$level="小心使用";
		}
		elsif($hash2{$cut[2]}<=1 and $cut[4] eq "增加用药"){
			$num2++;
			$level="增加剂量";
		}
		elsif($hash2{$cut[2]}<=1 and $cut[4] eq "小心使用"){
			$num3++;
			$level="小心使用";
		}
		if($hash2{$cut[2]}<=1){
			my @val=($cut[2],$level); 
			if(exists $hash{$cut[0]}){
				push @{$hash{$cut[0]}},[@val];
			}else{
				my @arr;
				$hash{$cut[0]}=\@arr;
				push @{$hash{$cut[0]}},[@val];
			}
		}
	}
	open IN,"<:encoding(utf-8)","$file1" || die "can't open the $file1!";    #华法林与氯吡格雷药物文件
	while(<IN>){
		chomp;
		my @cut=split/\t/;
		next if($_=~/^sample/);
		my $level;
		$hash2{$cut[9]}++;
		$cut[12]=~s/\s//g;
		$cut[8]=~s/心血管类/心血管病常用药/g;
		if($hash2{$cut[9]}<=1 and $cut[12] eq "正常使用"){
			$level="正常使用";
		}
		elsif($hash2{$cut[9]}<=1 and $cut[12] eq "建议换药"){
			$num1++;
			$level="不建议使用";
		}
		elsif($hash2{$cut[9]}<=1 and $cut[12] eq "调整剂量"){
			$num2++;
			$level="调整剂量";
		}
		if($hash2{$cut[9]}<=1){
			my @val=($cut[9],$level);
			if(exists $hash{$cut[8]}){
				push @{$hash{$cut[8]}},[@val];
			}else{
				my @arr;
				$hash{$cut[8]}=\@arr;
				push @{$hash{$cut[8]}},[@val];
			}
		}
	}
	my @med_class=("β受体阻断剂","利尿剂","血管紧张素转换酶抑制剂ACEI","血管紧张素II受体拮抗剂ARB","钙通道阻断剂","心血管病常用药","降脂药");
	close IN1;		
	close IN;
my $parthead;
if ($module1 eq 'Druga'){$parthead="检测结果";}
else{$parthead="心血管用药指导基因检测结果";}
my $biao=&mianpart($part1,$parthead);
my $head=<<'END';
\pagestyle{fancy}
\begin{YMpart}

END

my $title=<<'END';
\begin{center} \zihao{3}{\color{DarkBlue} \sye \sym 检测结果列表} \end{center} \\
\color{MyFontGray}
\arrayrulecolor{DarkBlue}
%\renewcommand\arraystretch{1.5}
\noindent                                                                                              
\fontsize{9.5}{13.3}\selectfont
END

my $info=qq[{此次药物敏感性分析一共为您检测了23种心血管病常用药，涵盖五大类降压药、降脂药和抗凝药。其中有$num1 种药物\\textcolor{red}{\\sym{不建议使用}}，有$num2 种药物建议您\\textcolor{orange}{\\sym{调整剂量}}，有$num3 种药物建议您\\textcolor{yellow}{\\sym{小心使用}}。具体情况参见下表和药物详细解析。}\\par];

my $table=<<'END';
\noindent
%\includegraphics[height=0.4cm]{13-2.pdf} \color{MyGreen}\fontsize{14}{14}{\sym{\ 检\,测\,结\,果\,列\,表}}\\
\begin{center}   
\fontsize{9}{9}
\color{MyFontGray}
\arrayrulecolor{DarkBlue}
\begin{VerbatimOut}{\jobname.vrb}                                                                                              
\renewcommand\arraystretch{1}
\begin{longtable}{ p{.32\textwidth} p{.36\textwidth} p{.25\textwidth} }
\rowcolor{DarkBlue}
\Gape[6pt]{\fontsize{12}{14.4}{\textcolor{white}{\sym{系统分类}}}} &\fontsize{12}{14.4}{\textcolor{white}{\sym{药物名称}}} & \fontsize{12}{14.4}{\textcolor{white}{\sym{建议等级}}} \tabularnewline 
END
my $tablejiewei=<<'END';
\end{longtable}
\end{VerbatimOut}
\LTXtable{\textwidth}{\jobname.vrb}
\end{center}
END


my $jiewei=<<'END';
\vspace*{-1cm}
\textcolor{green}{\sym{注：}}\\
\vspace*{-1.5cm}	
\begin{center}
\fontsize{9.5}{11.4}\selectfont
\color{MyFontGray}
\begin{tabular}{ p{.5\textwidth} p{.5\textwidth}}
\includegraphics[height=0.35cm]{药敏正常使用.pdf} 代表您能\textcolor{green}{\sym{正常使用}}该药物；& \includegraphics[height=0.35cm]{药敏小心使用.pdf} 代表该药物可能对您有副作用，建议\textcolor{yellow}{\sym{小心使用}}；\tabularnewline
\includegraphics[height=0.35cm]{药敏调整剂量.pdf} 代表药物剂量可能对您的影响较大，建议酌情\textcolor{orange}{\sym{调整剂量}}；& \includegraphics[height=0.35cm]{药敏建议换药.pdf} 代表该药物可能对您有严重不良反应，\textcolor{red}{\sym{不建议使用}}。\tabularnewline
\end{tabular}
\end{center}


%%%%%%%%%%%%%药物介绍 待加入

END

	my $tableinfo1=&tableinfor(@med_class);
	my $Latex="$biao"."$head"."$title"."$info"."$table"."$tableinfo1"."$tablejiewei"."$jiewei";
	return $Latex;
}
##################子程序################

sub tableinfor{
	my @med_class=@_;
	my $number=0;
	my ($color,$tinfo,$tinfo1,$tinfo2,$tableinfo);
	foreach my $a(@med_class){
if (exists $hash{$a}){
		my @c=@{$hash{$a}};
		my $number1=$number+1;
		my $n=@c;
		$number=$number+$n;
#		my $number1=$number+1;
#print "$n\n";
#print "$a\t$n\t$number\t$number1\n";		
		if($number==$#c+1){
			my $i;
			for($i=0;$i<=$#c;$i++){	
				my @arr=&table($a,$i,\@c);
				$tableinfo.="$arr[1]";
			}
		}
		if($number!=$#c+1 and $number1 % 2 ==1){
			my $i;
            for($i=0;$i<=$#c;$i++){ 
	           my @arr=&table($a,$i,\@c);
	           $tableinfo.="$arr[1]";
		    }
		}
		if($number1 % 2==0 and $number!=$#c+1){
			my $i;    
			for($i=0;$i<=$#c;$i++){
				my @arr=&table($a,$i,\@c);
				$c[$i]->[1]=~s/\s//g;          #有些建议等级后面有空格，解决‘药敏$c[$i]->[2].pdf’这个问题
				$color=$arr[0];
				my $ym=$arr[2];
				if($i==0 and $n!=1){
#print "$number1\t$a\t$c[$i]->[0]\n";
					$tinfo=qq[\\multirow{$n}{*}{\\thead{\\fontsize{9}{10.8}{$a}}} & \\rowcolor{LightBlue} $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline\n];
					$tableinfo.="$tinfo";
				}
				if($i==0 and $n==1){
					$tinfo=qq[$a & \\rowcolor{LightBlue} $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline \\hline\n];
                    $tableinfo.="$tinfo";
				}
				elsif(0<$i and $i<$#c and $i % 2 ==0){
					$tinfo1=qq[ & \\rowcolor{LightBlue} $c[$i]->[0]  & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline\n];
					$tableinfo.="$tinfo1";
				}
				elsif(0<$i and $i<$#c and $i % 2 ==1){
					$tinfo1=qq[ & $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline\n];
					$tableinfo.="$tinfo1";
				}
				elsif($i==$#c and $i % 2 ==0){
					$tinfo2=qq[ & \\rowcolor{LightBlue} $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline \\hline\n];
					$tableinfo.="$tinfo2";
				}
				elsif($i==$#c and $i % 2 ==1){
					$tinfo2=qq[ & $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline \\hline\n];
					$tableinfo.="$tinfo2";
				}
			}
		}
	}
 }                   #后加的
	return $tableinfo;
}


##########################子程序###########################
sub table{
	my ($a,$i,$b)=@_;
	my @c=@$b;
	my $n=$#c+1;
	my ($color,$tinfo,$tinfo1,$tinfo2,$tableinfo,$ym);
	if($c[$i]->[1] =~ /不建议使用/){   
		$color='red';
		$ym='建议换药';
	}
	elsif($c[$i]->[1] =~ /增加剂量/){
		$color='orange';
		$ym='调整剂量';
    }
	elsif($c[$i]->[1] =~ /调整剂量/){
		$color='orange';
		$ym='调整剂量';
	}
	elsif($c[$i]->[1] =~ /小心使用/){
		$color='yellow';
		$ym='小心使用';
	}
	elsif($c[$i]->[1] =~ /正常使用/){
		$color='green';
		$ym='正常使用';
	}
	$c[$i]->[1]=~s/\s//g;          #有些建议等级后面有空格，解决‘药敏$c[$i]->[2].pdf’这个问题
	if($i==0 and $n!=1){
#print "$c[$i]->[0]\n";
		$tinfo=qq[\\multirow{$n}{*}{\\thead{\\fontsize{9}{10.8}{$a}}} & $c[$i]->[0]  & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline\n];    
		$tableinfo.="$tinfo";
	}
	if($i==0 and $n==1){
		 $tinfo=qq[$a & $c[$i]->[0]  & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline \\hline \n];    
         $tableinfo.="$tinfo";
	}
	elsif(0<$i and $i<$#c and $i % 2 ==0){
		$tinfo1=qq[ & $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline\n];
		$tableinfo.="$tinfo1";
	}
	elsif(0<$i and $i<$#c and $i % 2 ==1){
		$tinfo1=qq[ & \\rowcolor{LightBlue} $c[$i]->[0]  & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline\n];
		$tableinfo.="$tinfo1";
	}
	elsif($i==$#c and $i % 2 ==0){
		$tinfo2=qq[ & $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline \\hline\n];
		$tableinfo.="$tinfo2";
	}
	elsif($i==$#c and $i % 2 ==1){
		$tinfo2=qq[ & \\rowcolor{LightBlue} $c[$i]->[0] & \\textcolor{$color}{\\sym{$c[$i]->[1]}} \\includegraphics[height=0.28cm]{药敏$ym.pdf} \\tabularnewline \\hline\n];
		$tableinfo.="$tinfo2";
	}
	my @array=($color,$tableinfo,$ym);
	return @array;
}
sub mianpart{
	my @input=@_;
	my $cover=
qq[%%%%%%%%%%%%%%%%% 心血管用药 结果 %%%%%%%%%%%%%%%

\\newpage
\\addcontentsline{toc}{section}{Part$input[0] — $input[1]}
\\addcontentsline{toc}{subsection}{ⅰ 检测结果}
\\ThisTileWallPaper{\\paperwidth}{\\paperheight}{XKhead.pdf}
\\begin{center} \\zihao{2}{\\color{DarkBlue} \\sye \\sym Part$input[0] — $input[1]}} \\end{center} \\\\];
return $cover;
}

}
1;
