use strict;
use warnings;
use utf8;
use Encode;
use Getopt::Long;
use POSIX;
use Data::Dumper;
use File::Basename;
use Smart::Comments;
use Getopt::Long;
use Spreadsheet::XLSX;
use List::Util qw/max min sum maxstr minstr shuffle/;

my ($samples,$gresult,$yresult,$dresult,$out,$usage)=('','','','','','');
GetOptions(
    'sample=s'=> \$samples,
    'gresult=s' => \$gresult,
    'yresult=s' => \$yresult,
    'dresult=s' => \$dresult,
    'out=s'  => \$out,
);

$usage=<<END;

perl $0
        -sample 样本ID
        -gresult,-yresult,-dresult   用药结果文件
        -out    样本报告输出路径
END

if($samples eq '' or $gresult eq '' or $yresult  eq '' or $dresult eq '')
{
        die "$usage\n";
}

$out||= "./";
my $bin = "/PUBLIC/work/panqi/diagnosis_report/script";
my $hutongbiao="/PUBLIC/pipline/database/sheet_hutong/";
my $drugclassfile="/PUBLIC/pipline/script/Report/XKFW/V2_171130/database/Drug.fenlei.txt";
my ($file,$file2) = judge_lastest_sitedb($hutongbiao);
`perl $bin/spectrum_XLSX.pl $file` unless(-e $file2);

my %durgresult;
my @drugpx;
my ($yemei,$biaoti,$info_sample,$info_samplet,$firstinfo,$num1,$num2,$num3,$drugtable1,$drugtable2,$parttail);
&read_drug_PX($drugclassfile);
&read_sample_info($file2);
&read_drug_result($gresult);
&read_drug_result($yresult);
&read_drug_result($dresult);
#&sub_table_summary;
&sub_table_xiangxi;
###################report tex#######
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
\newfontfamily\sye{SourceHanSansCN-Light.ttf}
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
\fancyhead[LO]{\includegraphics[height=1cm]{cardiologo.png}}
\fancyhead[LE]{\includegraphics[height=1cm]{cardiologo.png}}
\fancyfoot[CO,CE]{\color{MyFontGray} \fontsize{8}{9.6}\selectfont{北京诺禾心康基因科技有限公司 \ \ \ \ \ \ \ \ 地址：北京市海淀区金码大厦B
座16层 \ \ \ \ \ \ \ \  咨询电话：010-82837801-275}}
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
$biaoti="心血管用药指导基因检测";
$info_samplet=qq[\\section*{\\sym{$biaoti}}
\\zihao{-4} {\\sym{检测编号：$samples}} \\\\ \\\\
\\color{MyFontGray} \\fontsize{10.5}{12.6}\\selectfont
\\renewcommand\\arraystretch{1.5}
\\noindent
\\begin{tabular}{  L{3.5cm}  L{3.5cm}  L{3.5cm}  L{4cm}  } \\hline
];



$num1=$durgresult{number}{'不建议使用'};
$num2=$durgresult{number}{'调整剂量'}+$durgresult{number}{'增加剂量'};
$num3=$durgresult{number}{'小心使用'};
$firstinfo=qq[{此次药物敏感性分析一共为您检测了40种心血管病常用药，涵盖降压药、降脂药、降糖药和抗凝药。其中有$num1 种药物{\\sym{不建议使用}}，有$num2 种药物建议您{\\sym{调整剂量}}，有$num3 种药物建议您{\\sym{小心使用}}。具体情况参见下表和药物详细解析。}\\par];

&sub_table_summary;


$parttail=<<'END';
\zihao{6}
\vspace*{30mm}
检测机构：北京诺禾心康基因科技有限公司\\
公司地址：北京市海淀区金码大厦B座16层\\
联系电话：010-82837801-275\\
\end{document} 
END

my $texsample=$yemei.$info_samplet.$info_sample.$drugtable1.$drugtable2;
open(OUT2,">:utf8","$out/$samples.drug.tex");
print OUT2 $texsample;
close OUT2;
######################################################
sub read_sample_info{
	my ($file)=@_;
	open IN, "<:utf8","$file" or die $!;
	<IN>;
	while(<IN>)
	{
        	chomp;
	        my @cut=split /\t/, $_;
        	if($samples eq $cut[2])
        	{
                	my $name=$cut[4];
	                my $sex=$cut[7];
        	        my $year=$cut[8];
                	$cut[1]=~s/\//-/g;
	                my $id_len=length $cut[3];
			my ($shizi,$info1,$jiazu,$clinical_info)=("口腔拭子","","-","");
        	        if($id_len > 8){$shizi="静脉血";}
                	if($shizi eq "口腔拭子")
                	{
                        	$info1=qq[\\ \\ \\rlap{\${\\surd}\$}\${\\square}\$ 口腔拭子 \\ \\ \${\\square}\$ 静脉血];
                	}
	                elsif($shizi eq "静脉血")
        	        {
                	        $info1=qq[\\ \\ \${\\square}\$ 口腔拭子 \\ \\ \\rlap{\${\\surd}\$}\${\\square}\$ 静脉血];
                	}

                	$info_sample=qq[{\\sym{姓名：}} $name & {\\sym{性别：}} $sex & {\\sym{年龄：}} $year & {\\sym{收样日期：}} $cut[1] \\\\ \\hline
\\multicolumn{2}{l} {{\\sym{临床诊断：}} $clinical_info } & \\multicolumn{2}{l} {{\\sym{家族史：}} $jiazu } \\\\ \\hline
\\multicolumn{2}{l} {{\\sym{样本类型：}} $info1} & \\multicolumn{2}{l} {{\\sym{检测方法：}}高通量测序} \\\\ \\hline
\\end{tabular}
];

        	}
	}
	close IN;
}
#########################################################3
sub read_drug_PX{
	my ($file)=@_;
	open IN,"<:utf8","$file" || die "can't open the $file!";
	while(<IN>){
                chomp;
                my @cut=split/\t/,$_;
		push @drugpx,[$cut[2],$cut[1],$cut[0]];
	}
	close IN;
}

sub read_drug_result{
	my ($file)=@_;
        my @path=split/\//,$file;
        my $samplefile=(split/\./,$path[-1])[0];
        die "样本ID和用药文件ID不相符" if($samplefile ne $samples);
        my %head;
        my %drug_flag;
	open IN,"<:utf8","$file" || die "can't open the $file!"; 
	while(<IN>){
		chomp;
		my @cut=split/\t/,$_;
                if(/^药物类别/ || /^#/ || /^sample/i){
			for(my $i=0;$i<@cut;$i++){
				$head{$cut[$i]}=$i;
			}	
		}else{
			my $drugclassi=$head{'药物类别'};    my $drugclass=$cut[$drugclassi];
                        my $drugnamei=(exists $head{'中文名称'})?$head{'中文名称'}:$head{'药物名称'};     my $drugname=$cut[$drugnamei];
			my $jianyileveli=$head{'建议等级'};  my $jianyilevel=$cut[$jianyileveli]; 
                        my $rsnumberi=(exists $head{'Rs_Number'})?$head{'Rs_Number'}:$head{'Rs_number'};  my $rsnumber=$cut[$rsnumberi];
                        my $gnamei=(exists $head{'Gene'})?$head{'Gene'}:$head{'基因'};                    my $gname=$cut[$gnamei];
                        my $genotypei=(exists $head{'检测基因型'})?$head{'检测基因型'}:$head{'基因型'};       my $genotype=$cut[$genotypei];
                        my $siteresulti=$head{'位点结论'};   my $siteresult=$cut[$siteresulti];
                        my $othernamei=$head{'其他名称'};    my $othername=$cut[$othernamei];
                        my $drugsuggesti=$head{'用药建议'};  my $drugsuggest=$cut[$drugsuggesti];
                        #print "$samples\t$drugclass\t$drugname\t$jianyilevel\t$rsnumber\t$gname\t$genotype\t$siteresult\t$othername\t$drugsuggest\n";
		######################
			my $levelfinal='';
                        $jianyilevel=~s/\s+$//;
			if($jianyilevel eq "正常用药" || $jianyilevel eq "正常使用"){$levelfinal='正常用药';}
                	elsif($jianyilevel eq "不建议使用" || $jianyilevel eq "建议换药"){$levelfinal='不建议使用';}
                        elsif($jianyilevel eq "增加用药、小心使用" || $jianyilevel eq "小心使用"){$levelfinal="小心使用";}
                        elsif($jianyilevel eq "增加用药"){$levelfinal='增加剂量';}
			elsif($jianyilevel eq "调整剂量"){$levelfinal='调整剂量';}
			elsif($jianyilevel eq "酌情使用"){$levelfinal='小心使用';}
                        if(!exists $drug_flag{$drugname} &&  $levelfinal ne  ''){
				$durgresult{table1}{$drugname}{'level'}=$levelfinal;
                                $durgresult{table1}{$drugname}{'sugg'}=$drugsuggest;
                                $durgresult{number}{$levelfinal}++;
                                $drug_flag{$drugname}=1;  
                        }
                        push @{$durgresult{table2}{$drugname}},[$gname,$rsnumber,$genotype,$siteresult];	
		}
	}
	close IN;
}

sub sub_table_summary{
	$drugtable1=qq[
\\noindent \\\\ \\\\ \\\\
\\vspace*{6mm}
\\zihao{-4}{\\sym{检测结果列表：}} \\\\
\\renewcommand\\arraystretch{1.5}
\\fontsize{9.5}{10.6}\\selectfont{\\color{MyFontGray}}
$firstinfo
\\begin{longtable}{| C{15mm} | C{25mm} | C{20mm} | C{20mm} | L{70mm}|}  
\\hline
{\\sym{药物大类}}  &  {\\sym{药物分类}}  &  {\\sym{药物名称}}   & {\\sym{建议等级}} & \\makecell[c]{\\sym{用药建议}}  \\\\\\hline
];
	for(my $i=0;$i<@drugpx;$i++){
		my $dclass=$drugpx[$i][0];  my $dsubc=$drugpx[$i][1]; my $dname=$drugpx[$i][2];
		if(exists $durgresult{table1}{$dname}{'level'}){
			  $drugtable1.="$dclass  \&  $dsubc  \&  $dname \&  $durgresult{table1}{$dname}{'level'}  \& $durgresult{table1}{$dname}{'sugg'} \\\\\\hline"."\n";
		}
	}
	$drugtable1.=<<'END';
\end{longtable} 
END
}
sub sub_table_xiangxi{
	$drugtable2=qq[
	\\noindent \\\\ \\\\
\\zihao{-4} {\\sym{检测结果详细信息}} \\\\
\\renewcommand\\arraystretch{1.5}
\\fontsize{9.5}{10.6}
\\begin{longtable}{ | C{3cm}  |  C{3cm}  |  C{1.5cm}  |  C{2.5cm}  | C{1cm}|  L{4cm} | } 
\\hline			
{\\sym{药物种类}} & {\\sym{药物名称}} & {\\sym{基因名称}} & {\\sym{位点信息}}  & {\\sym{基因型}}  & \\makecell[c]{\\sym{结果说明}}   \\\\\\hline
];	
	for(my $i=0;$i<@drugpx;$i++){
		my $dclass=$drugpx[$i][0];  my $dsubc=$drugpx[$i][1]; my $dname=$drugpx[$i][2];
		if(exists $durgresult{table1}{$dname}{'level'} && $durgresult{table1}{$dname}{'level'} ne '正常用药'){
			for(my $i=0;$i<@{$durgresult{table2}{$dname}};$i++){
                                 print "$dname\t${$durgresult{table2}{$dname}}[$i][2]\t${$durgresult{table2}{$dname}}[$i][3]\n";
                                if(${$durgresult{table2}{$dname}}[$i][3] ne '-'){
					$drugtable2.="\\sye{$dsubc} \&  $dname \& \\sye{ ${$durgresult{table2}{$dname}}[$i][0] } \& \\sye{ ${$durgresult{table2}{$dname}}[$i][1] } \& \\sye{ ${$durgresult{table2}{$dname}}[$i][2] } \& \\sye{ ${$durgresult{table2}{$dname}}[$i][3] }  \\\\\\hline "."\n";
				}
			}
		}
	}
	$drugtable2.=<<'END';
\end{longtable} 
\end{document}
END
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
                if($date=~/-(\d+)\.xlsx/)
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

