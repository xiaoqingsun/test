#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use utf8;
use Encode;
use Getopt::Long;
use Data::Dumper;
use Smart::Comments;
use Text::Iconv;
use Spreadsheet::XLSX;
use List::Util qw/max min sum maxstr minstr shuffle/;

my ($gxyresult,$csresult,$dmresult,$wesresult,$gxy_flag,$cs_flag,$dm_flag,$wes_flag,$sanger_flag,$sampleID,$outdir);
GetOptions(
        'gxy_file=s'=> \$gxyresult,
        'cs_file=s' => \$csresult,
	'dm_file=s' => \$dmresult,
	'wes_file=s'  => \$wesresult,
        'gxy_flag'  => \$gxy_flag,
	'cs_flag'  => \$cs_flag,
	'dm_flag'  => \$dm_flag,
	'wes_flag'  => \$wes_flag,
	'sanger_flag' => \$sanger_flag, 
	'sample=s' => \$sampleID,
        'outdir=s' => \$outdir,
);

#die "没有输入文件 " if (!$gxyresult && !$csresult && !$dmresult && !$wesresult);
$outdir ||= "./";
my $bin = "/PUBLIC/work/panqi/diagnosis_report/script";
my $hutongbiao="/PUBLIC/pipline/database/sheet_hutong/";
my $drugclassfile="/PUBLIC/pipline/script/Report/XKFW/V2_171130/database/Drug.fenlei.txt";
my $hospital="/PUBLIC/pipline/script/Report/XKFW/V3_180824/script/product.config.txt";
my ($file,$file2) = judge_lastest_sitedb($hutongbiao);
`perl $bin/spectrum_XLSX.pl $file` unless(-e $file2);

my %sample_info;
my %sitetables=();
my %site_genotype=();
my %site_info=();
my %sample_sex=();
my @sitetable=();
my @sangerperson=();
my %genotype_hash=('野生纯合'=>0,'变异杂合'=>0,'变异纯合'=>0,'半合子'=>0);
my ($yiyuan,$header,$footer)=('','','');

my $fig="$outdir/fig/$sampleID.family.site.txt";
if(-e $fig)
{
	&read_sample_info_PCR_yiyuan($file2);
	&read_sample_sex($file2);
	&read_sample_info_PCR($fig);
	&read_site_table_PCR($fig);
}
else{
	&read_sample_info($file2);
	&read_site_table($gxyresult) if($gxy_flag && $gxy_flag==1);
	&read_site_table($csresult) if($cs_flag && $cs_flag==1);
	&read_site_table($dmresult) if($dm_flag && $dm_flag==1);
	&read_site_table($wesresult) if($wes_flag && $wes_flag==1);
}

########################tex;
my ($titlewen,$table,$shuoming,$conclu,$weidian,$endwen);

&read_header($hospital);
$titlewen=&texhead;
$table=&textable;
$shuoming=<<'END';
%\par \vspace*{2mm}
\\备注：参考 hg19 p37版本的参考基因组 \\ \\
%\vspace*{2cm}
END
$weidian=&outputsite;
$conclu=<<'END';
{\sym \sye{ \zihao{4} 3.检测结果说明：}} \\ \\
END
if(-e $fig){
$conclu.=&texconclu;
}
$endwen="\\end{document}";

my $texsample=$titlewen.$table.$shuoming.$conclu.$weidian.$endwen;
open(OUT2,">:utf8","$outdir/$sampleID\_JXsanger.tex");
print OUT2 $texsample;
close OUT2;

################################## tex sub
sub outputsite{
	my $sitet;
        $sitet.="\\newpage\n\\par\n \\vspace*{2mm}\n";
	my $n=0;
	my $len = @sitetable * @sangerperson; #wph add 190118
	foreach my $p(@sangerperson){
		if($p->[3] eq "先证者"){
			$sitet.="$p->[2] \\quad $p->[3]\\\\ \n";
		}else{
			$sitet.="$p->[2] \\quad 先证者$p->[3] \\\\ \n";
		}
		for(my $i=0;$i<@sitetable;$i++){
			$sitet.="检测位点： $sitetable[$i][3]  $sitetable[$i][4]，$sitetable[$i][5] \\\\";
			$fig="$outdir/fig/$p->[0]-$sitetable[$i][0]-$sitetable[$i][1].jpg";
			if(-e $fig)
			{
				if($n==$len-1 or $n%4==0 and $n!=$len){
					$sitet.=" \\\\ \\centerline{\\includegraphics[width=14cm,height=3.5cm]{$outdir/fig/$p->[0]-$sitetable[$i][0]-$sitetable[$i][1].jpg}}} \n";
				}else{
					$sitet.=" \\\\ \\centerline{\\includegraphics[width=14cm,height=3.5cm]{$outdir/fig/$p->[0]-$sitetable[$i][0]-$sitetable[$i][1].jpg}} \\\\ \\\\ \n";
				}
			}
			else{
				$sitet.="对照序列 \\\\ 样本序列\\\\ \n \\par\\vspace{30mm}\n";
			}
			$n++;
                        if($n%4==0 and $len!=$n){
                                $sitet.="\\newpage\n\\par\n \\vspace*{2mm}\n";
                        }
		}	
	}
	return $sitet;
}

sub texconclu{
	my $conclu;
	my $concluall;
	my $genotyp_count;
	my $geneall;
	my $mutation;
	my $unmutation;
	my $father_flag=0;
	my $mother_flag=0;
	my %count=();
	my @family;
	foreach my $key (keys %genotype_hash){
		if($genotype_hash{$key}){
			$genotyp_count.="$genotype_hash{$key}个$key，";
		} 
	}
	$genotyp_count=~s/变异杂合/杂合变异/g;
	$genotyp_count=~s/变异纯合/纯合变异/g;
	$genotyp_count=~s/半合子/半合子变异/g;
	$geneall=join("、",grep {++$count{$_} < 2} (values %sitetables) ); #去重后的基因
	$concluall.="先证者样本中检测出${genotyp_count}存在于${geneall}基因上，变异位点详情如上述。本报告为先证者$sample_info{'先证者'}[2]的家系验证报告，一代测序结果显示";
	my $site_count= keys %sitetables;
	@family= keys %sample_info;
	if(grep {$_ eq "母亲"} @family){ $mother_flag=1	}
	if(grep {$_ eq "父亲"} @family){ $father_flag=1 }
	print "位点数目$site_count\n";
	if($site_count==1){
		my $site=(keys %sitetables)[0];
		my @mutys;
		my @unmutys;
		(my $chr,my $start) = split /_/,$site;
		if(not $mother_flag and not $father_flag){
			$concluall.="，";
			$conclu="由于先证者的父母未验证，无法确定变异的遗传来源，";
		}elsif($mother_flag and $father_flag){
			if($site_genotype{'母亲'}{$chr}{$start} eq $site_genotype{'先证者'}{$chr}{$start} and $site_genotype{'父亲'}{$chr}{$start} ne $site_genotype{'先证者'}{$chr}{$start} or $site_genotype{'先证者'}{$chr}{$start} eq "半合子" and $site_genotype{'母亲'}{$chr}{$start} eq "变异杂合" ){
				$concluall.="，先证者的母亲$sample_info{'母亲'}[2]携带该变异，";
				$concluall.="父亲$sample_info{'父亲'}[2]未携带该变异，";
				$conclu="可见先证者的变异遗传自母亲，";
			}elsif($site_genotype{"父亲"}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start}){
				$concluall.="，先证者的父亲$sample_info{'父亲'}[2]携带该变异，";
				$concluall.="母亲$sample_info{'母亲'}[2]未携带该变异，";
				$conclu="可见先证者的变异遗传自父亲，";
			}else{
				$concluall.="父亲$sample_info{'父亲'}[2]、母亲$sample_info{'母亲'}[2]未携带该变异，";
				$conclu="可见先证者的变异为denovo变异。";
			}
		}elsif($mother_flag and not $father_flag){
			if($site_genotype{"母亲"}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start} or $site_genotype{'先证者'}{$chr}{$start} eq "半合子" and $site_genotype{'母亲'}{$chr}{$start} eq "变异杂合"){
			$concluall.="，先证者的母亲$sample_info{'母亲'}[2]携带该变异，";
			}else{
				$concluall.="，先证者的母亲$sample_info{'母亲'}[2]未携带该变异，";
			}
			$conclu="由于父亲未验证，无法确定变异的遗传来源，";
		}elsif(not $mother_flag and $father_flag){
			if($site_genotype{"父亲"}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start}){
				$concluall.="，先证者的父亲$sample_info{'父亲'}[2]携带该变异，";
			}else{
				$concluall.="，先证者的父亲$sample_info{'父亲'}[2]未携带该变异，";
			}
			 $conclu="由于母亲未验证，无法确定变异的遗传来源，";
		}
		foreach my $relation (@family){
			next if($relation eq "先证者" or $relation eq "母亲" or $relation eq "父亲");
			if($site_genotype{$relation}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start}){
				push @mutys,$relation.$sample_info{$relation}[2];
			}else{
				push @unmutys,$relation.$sample_info{$relation}[2];
			}
		}
		if(@mutys != 0){
			$mutation=join("、",@mutys)."携带该变异，"; 
			$mutation=~s/_1//g;
		}
		if(@unmutys != 0){
			$unmutation=join("、",@unmutys)."未携带该变异，";
			$unmutation=~s/_1//g;
		}
	$concluall.=$mutation.$unmutation.$conclu.qq[具体请结合临床分析，以下为本次检测样本的测序峰图：\\\\ \n];
	}else{
		$concluall.=qq[：\\\\ \n];
		for(my $i=0;$i<@sitetable;$i++){
			(my $chr,my $start) = ($sitetable[$i][0],$sitetable[$i][1]);
			$concluall.="$sitetable[$i][3]基因上的变异（$sitetable[${i}][4]，$sitetable[${i}][5]）";
			my @mutys;
			my @unmutys;
			if(not $mother_flag and not $father_flag){
				$concluall.="，";
				$conclu="由于先证者的父母未验证，无法确定变异的遗传来源；";
			}elsif($mother_flag and $father_flag){
				if($site_genotype{'母亲'}{$chr}{$start} eq $site_genotype{'先证者'}{$chr}{$start} and $site_genotype{'父亲'}{$chr}{$start} ne $site_genotype{'先证者'}{$chr}{$start} or $site_genotype{'先证者'}{$chr}{$start} eq "半合子" and $site_genotype{'母亲'}{$chr}{$start} eq "变异杂合"){
					$concluall.="，先证者的母亲$sample_info{'母亲'}[2]携带该变异，";
					$concluall.="父亲$sample_info{'父亲'}[2]未携带该变异，";
					$conclu="可见先证者的变异遗传自母亲；";
				}elsif($site_genotype{"父亲"}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start}){
					$concluall.="，先证者的父亲$sample_info{'父亲'}[2]携带该变异，";
					$concluall.="母亲$sample_info{'母亲'}[2]未携带该变异，";
					$conclu="可见先证者的变异遗传自父亲；";
				}else{
					$concluall.="先证者的父亲$sample_info{'父亲'}[2]、母亲$sample_info{'母亲'}[2]未携带该变异，";
					$conclu="可见先证者的变异为denovo变异；";
				}
			}elsif($mother_flag and not $father_flag){
				if($site_genotype{"母亲"}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start} or $site_genotype{'先证者'}{$chr}{$start} eq "半合子" and $site_genotype{'母亲'}{$chr}{$start} eq "变异杂合"){
				$concluall.="，先证者的母亲$sample_info{'母亲'}[2]携带该变异，";
				}else{
					$concluall.="，先证者的母亲$sample_info{'母亲'}[2]未携带该变异，";
				}
				$conclu="由于父亲未验证，无法确定变异的遗传来源；";
			}elsif(not $mother_flag and $father_flag){
				if($site_genotype{"父亲"}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start}){
					$concluall.="，先证者的父亲$sample_info{'父亲'}[2]携带该变异，";
				}else{
					$concluall.="，父亲$sample_info{'父亲'}[2]未携带该变异，";
				}
			 	$conclu="由于母亲未验证，无法确定变异的遗传来源；";
			}
			foreach my $relation (@family){
				next if($relation eq "先证者" or $relation eq "母亲" or $relation eq "父亲");
				if($site_genotype{$relation}{$chr}{$start} eq $site_genotype{"先证者"}{$chr}{$start}){
					push @mutys,$relation.$sample_info{$relation}[2];
				}else{
					push @unmutys,$relation.$sample_info{$relation}[2];
				}
			}
			if(@mutys != 0){
				$mutation=join("、",@mutys)."携带该变异，"; 
				$mutation=~s/_1//g;
			}
			if(@unmutys != 0){
				$unmutation=join("、",@unmutys)."未携带该变异，";
				$unmutation=~s/_1//g;
			}
			$concluall.=$mutation.$unmutation.$conclu.qq[\\\\ \n];
			$mutation="";
			$unmutation="";
		}
		$concluall.=qq[具体请结合临床分析，以下为本次检测样本的测序峰图：\\\\  \n];
	}
	
	print $concluall;
	return $concluall;
}

sub textable{
	my $tb2="";
	my $tb=<<'END';
\par \vspace*{2mm}
\begin{center}
END
	$tb.=qq[\\zihao{2}{$sampleID\号家系一代测序结果} \\\\\n];

	$tb.=<<'END';
\end{center}
\par \vspace*{6mm}
\renewcommand\arraystretch{1.4}
\fontsize{12}{14.5}\selectfont 
{\sym \sye{ \zihao{4} 1.样本信息}} \\
\begin{longtable}{| C{40mm} | C{50mm} | C{35mm} | C{35mm} | } \hline
{\sym{姓名}} & {\sym{样品编号}} & {\sym{家庭关系}} & {\sym{样本类型}} \\ \hline
END
	foreach my $p(@sangerperson){
		$tb.=qq[$p->[2] & $p->[1] & $p->[3] & $p->[4] \\\\ \\hline \n];
	}
	
	$tb.=<<'END';
\end{longtable}
\vspace*{5mm}
{\sym \sye{ \zihao{4} 2.位点基因型的一代验证结果}} \\
%\tablefirsthead{
%\hline 样品编号  &  样品类型  &  检测位点   & 所在基因  & 检测结果  \\\hline }
%\tablehead{}
%\tabletail{\hline}
%\tablelasttail{}
%\begin{supertabular}{| m{20mm}<{\centering} | m{20mm}<{\centering} | m{55mm}<{\centering} | m{20mm}<{\centering} |m{30mm}<{\centering} | }
END

	my $num=0;
	my $num_4=0;
	my %num_sub_sample;
	my $nsite=@sangerperson;
	my $cycle;
	if( int($nsite/4)*4 == $nsite ){ $cycle	= int( $nsite/4 ); }
	if( int($nsite/4)*4 <  $nsite ){ $cycle = int( $nsite/4 ) + 1; }
	for(my $i=0;$i<$cycle;$i++)
	{
		#$tb2.="\\begin{longtable}{| C{55mm} | C{20mm} | ";
		if( ($i+1)*4 <= $nsite || $i==0)
		{
			$tb2.="\\begin{longtable}{| C{55mm} | C{20mm} | ";
			$tb2.="C{20mm} |" x 4;
		}
		else{
			$tb2.="\\begin{tabular}{| C{55mm} | C{20mm} | ";
			$tb2.="C{20mm} |" x 4;
		}
		$tb2.=" } \\hline\n";

		$tb2.="{\\sym{检测位点}}  & {\\sym{所在基因}} ";
		
		$num_4=0;
		while( $num<$nsite && $num<($i+1)*4 )
		{
			$tb2.="& {\\sym{\\makecell*[c]{$sangerperson[$num][2] \\\\ $sangerperson[$num][3]} }} ";
			$num_sub_sample{$num_4}=$sangerperson[$num][0];
			$num_4++;
			$num++;
		}
		$tb2.=" & " x (($i+1)*4-$num);
		$tb2.=" \\\\ \\hline \n";

		for(my $j=0;$j<@sitetable;$j++){
			my $site="chr$sitetable[$j][0]:$sitetable[$j][1]";
			if($sitetable[$j][2]=~/^rs/){ $site=$sitetable[$j][2]; }
			$tb2.="\\makecell*[c]{$site\\\\\($sitetable[$j][4],$sitetable[$j][5]\)}  \& $sitetable[$j][3] ";
			
			if(-e $fig)
			{
				if( ($i+1)*4 <= $nsite ){
					$tb2.="\& \\makecell*[c]{ $site_genotype{ $num_sub_sample{0} }{ $sitetable[$j][0] }{ $sitetable[$j][1]} }  \& \\makecell*[c]{ $site_genotype{ $num_sub_sample{1} }{ $sitetable[$j][0] }{ $sitetable[$j][1]} }  \& \\makecell*[c]{ $site_genotype{ $num_sub_sample{2} }{ $sitetable[$j][0] }{ $sitetable[$j][1]} }  \& \\makecell*[c]{ $site_genotype{ $num_sub_sample{3} }{ $sitetable[$j][0] }{ $sitetable[$j][1]} }";
				}
				else{
					for(my $k=0;$k<$num_4;$k++)
					{
						$tb2.="\& \\makecell*[c]{ $site_genotype{ $num_sub_sample{$k} }{ $sitetable[$j][0] }{ $sitetable[$j][1]} } ";
					}
					$tb2.=" \&"x(($i+1)*4-$nsite);
				}				
			}
			else{
				if( ($i+1)*4 <= $nsite ){$tb2.="\& \\makecell*[c]{- \\hspace*{1mm} - \\hspace*{1mm}- \\hspace*{1mm} -} " x 4;}
				else{$tb2.="\& \\makecell*[c]{- \\hspace*{1mm} - \\hspace*{1mm}- \\hspace*{1mm} -} " x ( $nsite-$i*4 )." \&" x(($i+1)*4-$nsite);}
			}
			
			$tb2.=" \\\\ \\hline \n";
		}
		$tb2.="\\end{longtable} \n \\vspace*{2mm}\n" if( ($i+1)*4 <= $nsite || $i==0);
		$tb2.="\\end{tabular}   \n \\vspace*{2mm}\n" if( ($i+1)*4 >  $nsite && $i!=0);
	}
	$tb.=$tb2;
	return $tb;

=head
	foreach my $p(@sangerperson){
		#if($p->[1]=~/NKH/){ $p->[1]= substr($p->[1],0,8)."\\\\".substr($p->[1],8);}#wph add 180929
		print "$p->[0]\t$p->[1]\t$p->[2]\t$p->[3]\n";
		if($nsite==1){
			my $site="chr$sitetable[0][0]:$sitetable[0][1]";
                        if($sitetable[0][2]=~/^rs/){
				$site=$sitetable[0][2];
			}
			$tb.="\\makecell*[c]{$p->[2]\\\\\($p->[1]\)\\\\$p->[3]} \& $p->[4] \&  \\makecell*[c]{$site\\\\\($sitetable[0][4],$sitetable[0][5]\)}  \& $sitetable[0][3]  \&  \\\\ \\hline\n";
		}else{
			for(my $i=0;$i<@sitetable;$i++){
                                my $site="chr$sitetable[$i][0]:$sitetable[$i][1]";
				if($sitetable[$i][2]=~/^rs/){
					$site=$sitetable[$i][2];
				}
				if($i==0){
					$tb.="\\multirow{$nsite}{*}{\\makecell*[c]{$p->[2]\\\\\($p->[1]\)\\\\$p->[3]}} \& \\multirow{$nsite}{*}{$p->[4]} \&  \\makecell*[c]{$site\\\\\($sitetable[$i][4],$sitetable[$i][5]\)}  \& $sitetable[$i][3]  \&  \\\\ \\cline{3-5}\n";
				}elsif($i==$nsite-1){
					$tb.="\& \& \\makecell*[c]{$site\\\\\($sitetable[$i][4],$sitetable[$i][5]\)}  \& $sitetable[$i][3]  \&  \\\\ \\hline\n";
				}else{
					$tb.="\& \& \\makecell*[c]{$site\\\\\($sitetable[$i][4],$sitetable[$i][5]\)}  \& $sitetable[$i][3]  \&  \\\\ \\cline{3-5}\n";
				}
			}
		}
	}
	$tb.="\\end{supertabular}\n";
=cut

}

sub texhead{
	my $title=<<'END';
%-*- coding: UTF-8 -*-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[UTF8,noindent,twoside]{ctexrep}
\usepackage{blkarray,fancyvrb,ltxtable,multirow,colortbl,array,makecell,fancyhdr,graphicx,multicol,balance,amssymb,setspace,dcolumn,booktabs,tabularx,longtable,wallpaper,tabu,makeidx,algorithm,algorithmic,changepage,fontspec,wrapfig,picinpar,enumerate,pdfpages,textcomp,supertabular}
\usepackage{flushend}
\usepackage[abs]{overpic}
\setlength{\parindent}{0pt} 
\usepackage[a4paper,left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}
\usepackage[table]{xcolor}
\usepackage[center]{titlesec}
\graphicspath{{/PUBLIC/work/panqi/diagnosis_report/Latex/figure/}}
%define colour
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
%Font setting
\setmainfont{SourceHanSansCN-Light.ttf}
\setsansfont{SourceHanSansCN-Light.ttf}
\setmonofont{SourceHanSansCN-Light.ttf}
\setCJKmainfont{SourceHanSansCN-Light.ttf}
\setCJKsansfont{SourceHanSansCN-Light.ttf}
\setCJKmonofont{SourceHanSansCN-Light.ttf}
\newfontfamily\sye{SourceHanSansCN-Medium.ttf}
\newCJKfontfamily\sym{SourceHanSansCN-Medium.ttf}
\newcolumntype{C}[1]{>{\centering\arraybackslash}m{#1}}
\newcolumntype{L}[1]{>{\raggedright\arraybackslash}m{#1}}
\newcolumntype{R}[1]{>{\raggedleft\arraybackslash}m{#1}}
\newcolumntype{Y}{>{\centering\arraybackslash}X}
\newcommand{\tabincell}[2]{\begin{tabular}{@{}#1@{}}#2\end{tabular}}
\titleformat{\section}{\centering\zihao{-2}}{}{2em}{}
\newcommand\mgape[1]{\gape{$\vcenter{\hbox{#1}}$}}  %%%表格中图片居中
\pagestyle{fancy}
\fancyhf{}
\renewcommand\headrulewidth{0.6pt}
\renewcommand{\footrulewidth}{0.0pt}
END
$title.=qq[$header];
$title.=qq[$footer];
$title.=<<'END';
%\fancyhead[LO]{\includegraphics[height=1cm]{cardiologo.png}}
%\fancyhead[LE]{\includegraphics[height=1cm]{cardiologo.png}}
%\fancyfoot[CO,CE]{\color{MyFontGray} \fontsize{8}{9.6}\selectfont{北京诺禾心康基因科技有限公司 \ \ \ \ \ \ \ \ 地址：北京市朝阳区酒仙桥电子城IT产业园C区C1楼A座6层 \ \ \ \ \ \ \ \  咨询电话：4000-801-808}}
\makeatletter %双线页眉
\def\headrule{{\color{DarkBlue}\if@fancyplain\let\headrulewidth\plainheadrulewidth\fi
\hrule\@height 0.3pt \@width\headwidth\vskip1pt%上面线为1pt粗
\hrule\@height 3.0pt\@width\headwidth %下面0.5pt粗
\vskip-2\headrulewidth\vskip-3pt} %两条线的距离1pt
\vspace{3mm}} %双线与下面正文之间的垂直间距
\begin{document}%%%%%%%%%%%%%%%%%%%%%%% Report Face %%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage
\vspace*{-4mm}
END
	return $title;
}
#############################
sub read_sample_info_PCR_yiyuan{
	my $file=shift;
	my $id=$sampleID;
	open(IN,"<:utf8",$file) ||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		if($t[2] eq $sampleID){
			$yiyuan=$t[10];
		}
	}
	close IN;
}

sub read_sample_info_PCR{
	my $file=shift;
	my %hash;
	open(IN,"<:utf8",$file) ||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		next if($t[1] ne $sampleID);
		if(!exists $hash{$t[2]})
		{
			$t[0]=~s/妈妈/母亲/g;
			$t[0]=~s/爸爸/父亲/g;
			push @sangerperson,[$t[2],$t[11],$t[10],$t[0],$t[12]];
			if(!exists $sample_info{$t[0]}){
				$sample_info{$t[0]}=[$t[2],$t[11],$t[10],$t[0],$t[12]]; #编号，拭子，姓名，关系，样本类型
			}else{
				$sample_info{"$t[0]_1"}=[$t[2],$t[11],$t[10],$t[0],$t[12]];
			}
			$hash{$t[2]}=1;
		}
	}
	close IN;
}

sub read_sample_sex{ #wph add 190508 添加性别
	my ($file)=@_;
	my $id=$sampleID;
	$id=~s/^N(0+)//;
	open(IN,"<:utf8",$file) ||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		if($t[2] eq $sampleID){
			$sample_sex{'先证者'} = $t[7];
		}else{
			if($t[6]=~/一代验证-(\d+)-(.*)/){
				if($1 eq $id){
					my $guanxi=$2;
					$guanxi=~s/妈妈/母亲/g;
					$guanxi=~s/爸爸/父亲/g;
					if(!exists $sample_sex{$guanxi}){
						$sample_sex{$guanxi} = $t[7];
					}else{
						$sample_sex{"$guanxi\_1"} = $t[7];
					}
				}
			}
		}
	}
}

sub read_sample_info{
	my ($file)=@_;
	my $id=$sampleID;	
	$id=~s/^N(0+)//;
	open(IN,"<:utf8",$file) ||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		my $type=$t[11];
		my $id_shizi=(split /\//,$t[3])[0];
		if($id_shizi=~/([0-9]+)/)            ### changed by panqi 20181108 
		{
			if(length($1) == 8 && $type=~/0|None/)
			{
				$type="口腔拭子";
			}
		}
		if($t[2] eq $sampleID){
			$yiyuan=$t[10];
			#push @sangerperson,[$t[2],$t[3],$t[4],'先证者',$t[11]];
			push @sangerperson,[$t[2],$t[3],$t[4],'先证者',$type];
			$sample_info{'先证者'} = [$t[2],$t[3],$t[4],'先证者',$type];
			$sample_sex{'先证者'} = $t[7];
		}else{
			if($t[6]=~/一代验证-(\d+)-(.*)/){
				if($1 eq $id){
					#push @sangerperson,[$t[2],$t[3],$t[4],$2,$t[11]];
					$2=~s/妈妈/母亲/g;
					$2=~s/爸爸/父亲/g;
					push @sangerperson,[$t[2],$t[3],$t[4],$2,$type];
					$sample_info{$2} = [$t[2],$t[3],$t[4],$2,$type];
					$sample_sex{$2} = $t[7];
				}
			}
		}
	}
	close IN;
}

sub read_header{
	my ($file)=@_;
	open(IN,"<:utf8",$file) ||die;
	while(<IN>){
		next if(/^product/);
		chomp;
		my @cut=split /\t/,$_;
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
	close IN;

}

sub read_site_table{
	my ($file)=@_;
	open(IN,"<:utf8",$file) ||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		if($t[0] eq $sampleID && $t[1] ne '-'  &&  !exists $sitetables{"$t[1]_$t[2]"}){ #2,3,15,86-88
			$sitetables{"$t[1]_$t[2]"} = $t[85]; #赋值为基因名字
			if($t[87]=~/_/ or $t[86]=~/_/){$t[86]=~s/_/\\_/;$t[87]=~s/_/\\_/;  #wph changed 181107
			print "###$t[86]\t$t[87]\n";}                                      #wph changed 181107
                        print "$t[1],$t[2],$t[18],$t[85],$t[86],$t[87]\n";                 ### panqi 180817
			push @sitetable,[$t[1],$t[2],$t[18],$t[85],$t[86],$t[87]];	   ##chr,start,rsnum,gene,cds,pep
		}
	}
	close IN;
}

sub read_site_table_PCR{
	my $file=shift;
	my %het_hom = ('00'=>'野生纯合','01'=>'变异杂合','11'=>'变异纯合');
	open(IN,"<:utf8",$file) ||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		$t[0]=~s/妈妈/母亲/g;
		$t[0]=~s/爸爸/父亲/g;
		if($t[2] eq $sampleID && $t[2] ne '-'  &&  !exists $sitetables{"$t[3]_$t[4]"}){
			my $site_type = $het_hom{$t[15]};
			if($t[3] eq "X" and $sample_sex{$t[0]} eq "男" and $site_type eq "变异纯合"){
				$site_type = "半合子";
			}
			$genotype_hash{$site_type} ++;
			$sitetables{"$t[3]_$t[4]"} = $t[6]; #赋值为基因名字
			if($t[8]=~/_/ or $t[7]=~/_/){$t[8]=~s/_/\\_/;$t[7]=~s/_/\\_/;print "###$t[7]\t$t[8]\n";}
			print "$t[3],$t[4],$t[5],$t[6],$t[7],$t[8]\n";
			push @sitetable,[$t[3],$t[4],$t[5],$t[6],$t[7],$t[8]];		##chr,start,rsnum,gene,cds,pep
		}
		my $site_type = $het_hom{$t[15]};
		if($t[3] eq "X" and $sample_sex{$t[0]} eq "男" and $site_type eq "变异纯合"){
			$site_type = "半合子";
		}
		$site_genotype{$t[2]}{$t[3]}{$t[4]}=$site_type;
		if(exists $site_genotype{$t[0]}{$t[3]}{$t[4]}){
			$site_genotype{"$t[0]_1"}{$t[3]}{$t[4]}=$site_type;
		}else{
			$site_genotype{$t[0]}{$t[3]}{$t[4]}=$site_type;
		}
	}
	close IN;
}

sub judge_lastest_sitedb{
        my $product_site=shift;
        my @all_date;
        my $lastest_date;
        my $date;

        my @all_path=glob "$product_site/2019心*";
        foreach(@all_path)
        {
                $date=(split /\//,$_)[-1];
                if($date=~/-(\d+)\.xlsx/)
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

