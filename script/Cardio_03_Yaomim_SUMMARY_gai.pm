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
use encoding "utf-8";
#use utf8;
use Encode;
use POSIX;
use Data::Dumper;


my %hash;
my %yhash;
my %DaYaclass;

sub summary{
	my ($gresult,$yresult,$dresult,$part,$module)=@_;
############################### main ##############################
    my $tex=&shuru($gresult,$yresult,$dresult,$part,$module);
	return $tex;

#################读入文件进行处理##############
sub shuru{
	my $file=shift;
	my $file1=shift;
	my $file2=shift;
	my $part1=shift;
    	my $module1=shift;
	my @path=split/\//,$file;
	my @sample_path=split/\./,$path[-1];
	my $sample=$sample_path[0];
	my $num=0;my $num1=0;my $num2=0;my $num3=0;
	my (%hash1,%hash2,%hash3);
	my (%ghash,@detail);
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

		@detail=($cut[7],$cut[8],$cut[9],$cut[27],$cut[30]); #基因名 位点 基因型 临床注释 用药建议
		push (@{$ghash{$cut[2]}},@detail);  #key:药物名称  value:[[基因名,位点,基因型,临床注释,用药建议]]

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
		if($cut[8]=~/心血管类/){ $cut[8]="心血管病常用药"; }
		#$cut[8]=~s/心血管类/心血管病常用药/g;
		@detail=($cut[6],$cut[7],$cut[11],"-",$cut[14]); #基因名 位点 基因型 临床注释 用药建议
		push (@{$ghash{$cut[9]}},@detail); #key:药物名称  value:[[基因名,位点,基因型,临床注释,用药建议]]
		if($hash2{$cut[9]}<=1 and $cut[12] eq "正常使用"){
			$level="正常使用";
		}
		elsif($hash2{$cut[9]}<=1 and $cut[12] eq "建议换药"){
			$num1++;
			$level="不建议使用";
		}
		elsif($hash2{$cut[9]}<=1 and $cut[12] eq "*17:0/1;建议换药"){   ### panqi add 180924
			$num1++;
			$level="不建议使用;*17:0/1";
		}
		elsif($hash2{$cut[9]}<=1 and $cut[12] eq "*17:1/1;建议换药"){
			$num1++;
			$level="不建议使用;*17:1/1";
		}
		elsif($hash2{$cut[9]}<=1 and $cut[12] eq "调整剂量"){
			$num2++;
			$level="调整剂量";
		}
		elsif($hash2{$cut[9]}<=1 and $cut[12] eq "酌情使用"){
			$num3++;
			$level="小心使用";
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

  open IN2,"<:encoding(utf-8)","$file2" || die "can't open the $file1!";    #华法林与氯吡格雷药物文件
	while(<IN2>){
		chomp;
		my @cut=split/\t/;
		next if($_=~/^sample/);
		my $level;
		$hash2{$cut[9]}++;
		$cut[12]=~s/\s//g;
		$cut[8]=~s/心血管类/心血管病常用药/g;
		@detail=($cut[6],$cut[7],$cut[11],"-",$cut[14]); #基因名 位点 基因型 临床注释 用药建议
		push (@{$ghash{$cut[9]}},@detail); #key:药物名称  value:[[基因名,位点,基因型,临床注释,用药建议]]
		if($hash2{$cut[9]}<=1 and $cut[12] eq "正常使用"){
			$level="正常使用";
		}elsif($hash2{$cut[9]}<=1 and $cut[12] eq "建议换药"){
			$num1++;
			$level="不建议使用";
		}elsif($hash2{$cut[9]}<=1 and $cut[12] eq "调整剂量"){
			$num2++;
			$level="调整剂量";
		}elsif($hash2{$cut[9]}<=1 and $cut[12] eq "酌情使用"){
			$num3++;
			$level="小心使用";
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

	my @med_class=("利尿剂","β受体阻断剂","血管紧张素转换酶抑制剂ACEI","血管紧张素II受体拮抗剂ARB","钙通道阻断剂","心血管病常用抗凝药","双胍类药物","磷酸脲类促胰岛素分泌药","胰岛素增敏类药物","肠促胰岛素类似物","他汀类药物","贝特类药物","胆固醇吸收抑制剂");
  close IN2;
	close IN1;
	close IN;

my %Dalei;
my %Yalei;
my %alldrug;
open IN3,"<:encoding(utf-8)","/PUBLIC/pipline/script/Report/XKFW/V2_171130/database/Drug.fenlei.txt" || die "can't open the Drug.fenlei.txt!";  #大类 亚类 药物名称
        while(<IN3>){
                chomp;
                my @cut=split/\t/;
                next if($_=~/^药物名称/);
                $Dalei{$cut[0]}=$cut[2];  #key:药物名称 value:大类
		$Yalei{$cut[0]}=$cut[1];  #key:药物名称 value:亚类
}
close IN3;


foreach my $key1 (keys %hash){
      my @druglevel=@{$hash{$key1}};
		for(my $i=0;$i<scalar(@druglevel);$i++){
			$yhash{$Dalei{$druglevel[$i][0]}}{$Yalei{$druglevel[$i][0]}}{$druglevel[$i][0]} = $druglevel[$i][1]; #%yhash  key:药物名称  value:建议等级
			$alldrug{$druglevel[$i][0]}=1;
		}
}

open IN3,"<:encoding(utf-8)","/PUBLIC/pipline/script/Report/XKFW/V2_171130/database/Drug.fenlei.txt" || die "can't open the Drug.fenlei.txt!";  #大类 亚类 药物名称
 while(<IN3>){
                chomp;
                my @cut=split/\t/;
                next if($_=~/^药物名称/);
		next if(!exists $alldrug{$cut[0]});
                if(!exists $DaYaclass{$cut[2]}{ya_f}{$cut[1]}){
                        push @{$DaYaclass{$cut[2]}{ya_c}},$cut[1];
                         $DaYaclass{$cut[2]}{ya_f}{$cut[1]}=1;
                }
                push @{$DaYaclass{$cut[2]}{drug_c}{$cut[1]}},$cut[0];
}
close IN3;

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

my $info=qq[{此次药物敏感性分析一共为您检测了40种心血管病常用药，涵盖降压药、降脂药、降糖药和抗凝药。其中有$num1 种药物\\textcolor{red}{\\sym{不建议使用}}，有$num2 种药物建议您\\textcolor{orange}{\\sym{调整剂量}}，有$num3 种药物建议您\\textcolor{yellow}{\\sym{小心使用}}。具体情况参见下表和药物详细解析。}\\par];

my $table=<<'END';
\begin{spacing}{1.3}
\zihao{-5}
\color{MyFontGray}
\arrayrulecolor{DarkBlue}
\begin{VerbatimOut}{\jobname.vrb}
\renewcommand\arraystretch{1}
\begin{longtable}{L{.1\textwidth} L{.28\textwidth} L{.33\textwidth} L{.2\textwidth} }
\rowcolor{DarkBlue}
& \Gape[6pt]{\fontsize{12}{14.4}{\textcolor{white}{\sym{系统分类}}}} &\fontsize{12}{14.4}{\textcolor{white}{\sym{药物名称}}} & \fontsize{12}{14.4}{\textcolor{white}{\sym{建议等级}}} \tabularnewline
END
my $tablejiewei=<<'END';
\end{longtable}
\end{VerbatimOut}
\LTXtable{\textwidth}{\jobname.vrb}
\end{spacing}
END


my $jiewei=<<'END';
\vspace*{-0.5cm}
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

	my $tableinfo1=&tableinfor();
	#my $tableinfo2=&tableinfor2(%ghash);
	my $Latex="$biao"."$head"."$title"."$info"."$table"."$tableinfo1"."$tablejiewei"."$jiewei";
	return $Latex;
}
##################子程序################

sub tableinfor{
	my @alldalei=("高血压","2型糖尿病","高血脂","其他");
	my %dengji=(
        "正常使用"=>"\\textcolor{green}{\\sym{正常使用}} \\includegraphics[height=0.28cm]{药敏正常使用.pdf}",
        "增加剂量"=>"\\textcolor{orange}{\\sym{增加剂量}} \\includegraphics[height=0.28cm]{药敏调整剂量.pdf}",
        "调整剂量"=>"\\textcolor{orange}{\\sym{调整剂量}} \\includegraphics[height=0.28cm]{药敏调整剂量.pdf}",
        "不建议使用"=>"\\textcolor{red}{\\sym{不建议使用}} \\includegraphics[height=0.28cm]{药敏建议换药.pdf} ",
	"不建议使用;*17:0/1"=>"\\textcolor{red}{\\sym{不建议使用;*17:0/1}} \\includegraphics[height=0.28cm]{药敏建议换药.pdf} ",  ### panqi add 180924
	"不建议使用;*17:1/1"=>"\\textcolor{red}{\\sym{不建议使用;*17:1/1}} \\includegraphics[height=0.28cm]{药敏建议换药.pdf} ",
        "小心使用"=>"\\textcolor{yellow}{\\sym{小心使用}} \\includegraphics[height=0.28cm]{药敏小心使用.pdf}",
        );
	my $table="";
	my $num=0;
        foreach my $everydalei(@alldalei){
		my $yalei_num=0;
		my $multidalei=0;
		foreach my $ykey (@{$DaYaclass{$everydalei}{ya_c}}){
		#foreach my $ykey (keys %{$yhash{$everydalei}}){
			#my $yvalue=$yhash{$everydalei}{$ykey};
			my $drug_num=0;
			$num++;
			$yalei_num++;
			my $mulyalei=@{$DaYaclass{$everydalei}{drug_c}{$ykey}};
			foreach my $y2key (@{$DaYaclass{$everydalei}{drug_c}{$ykey}}){
			#foreach my $y2key (keys %{$yhash{$everydalei}{$ykey}}){
				my $y2value=$yhash{$everydalei}{$ykey}{$y2key};
				my $yalei;
				my $dengjitex;
				$drug_num++;
				$multidalei++;
				if (exists $dengji{$y2value}){
                        		$dengjitex=$dengji{$y2value};
                		}else{
                        		print "please check key $y2value is right";
               			}
				#my $mulyalei=keys%$yvalue;
				if($drug_num==$mulyalei && $yalei_num!=keys%{$yhash{$everydalei}}){
					$yalei=&get_rowcolor($num,$ykey,$mulyalei);
                        		$table.=" & $yalei & $y2key & $dengjitex \\tabularnewline"."\n";
				}elsif($drug_num==$mulyalei && $yalei_num==keys%{$yhash{$everydalei}}){
					$yalei=&get_rowcolor($num,$ykey,$mulyalei);
                        		$table.="\\multirow{-tihuan}{*}{\\zihao{5}\\sym{$everydalei}} & $yalei & $y2key & $dengjitex \\tabularnewline \\hline"."\n";
				}else{
					$yalei=&get_rowcolor($num,"",$mulyalei);
                                        $table.=" & $yalei & $y2key & $dengjitex \\tabularnewline"."\n";
				}
			}
		}
		$table=~s/tihuan/$multidalei/g;
	}
	return $table;
}

sub get_rowcolor{
	my $num=shift;
	my $yname=shift;
	my $multi=shift;
	my @arr;
	my $color="";
	my $yalei="";
	if($num%2==0){$color="\\rowcolor{LightBlue}"}
	if($yname eq ""){
		$yalei=$color;
	}else{
		if($multi==1){
			$yalei=$color.$yname;
		}else{
			$yalei=$color."\\multirow{-$multi}{3cm}{$yname}";
		}
	}
	return $yalei;
}

sub get_dengji{
	my $yname=shift;
	my $dengjitex;
	my %jianyi=(
	"正常使用"=>"\\textcolor{green}{\\sym{正常使用}} \\includegraphics[height=0.28cm]{药敏正常使用.pdf}",
	"增加剂量"=>"\\textcolor{orange}{\\sym{增加剂量}} \\includegraphics[height=0.28cm]{药敏调整剂量.pdf}",
	"调整剂量"=>"\\textcolor{orange}{\\sym{调整剂量}} \\includegraphics[height=0.28cm]{药敏调整剂量.pdf}",
	"不建议使用"=>"\\textcolor{red}{\\sym{不建议使用}} \\includegraphics[height=0.28cm]{药敏建议换药.pdf} ",
	"小心使用"=>"\\textcolor{yellow}{\\sym{小心使用}} \\includegraphics[height=0.28cm]{药敏小心使用.pdf}",
	);
	#while(my ($key,$value) = each(%yhash)){
	#	print "$key => $value\n";}
	if (exists $yhash{$yname}){
		my $level=$yhash{$yname};
		#$level=s/(^s+|s+$)//g;
		if (exists $jianyi{$level}){
			$dengjitex=$jianyi{$level};
		}else{
			print "please check key$level is right";
		}
	}else{
		print "please check key$yname is right";
	}
	return $dengjitex;
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
