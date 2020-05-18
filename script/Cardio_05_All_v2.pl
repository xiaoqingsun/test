#!/usr/bin/env perl
#===============================================================================
#
#         FILE: gxy_report.pl
#
#        USAGE: ./gxy_report.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (zhouwei),
# ORGANIZATION:
#      VERSION: 1.0
#      CREATED: 2017年04月06日 12时24分17秒
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
#use utf8;
use Encode;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use feature qw/say/;

my ($qc,$cfg,$in_dir,$right,$DM,$IDT,$gresult,$yresult,$dresult,$outdir,$cusi,$nothing,$bin,$type_new);
$in_dir='';
$right='';$cusi='';$gresult='';$yresult='';$dresult='';$nothing='';$bin='';

GetOptions(
    'cfg=s' => \$cfg,
    'in=s' => \$in_dir,
    'right:s'=> \$right,
    'cusi:s'=> \$cusi,  #${ID}_CS.pdf.txt
	'gresult:s'=> \$gresult,
	'yresult:s'=> \$yresult,
	'dresult:s'=> \$dresult,
	'module=s' => \$nothing,  #已无用
	'type=s'   => \$type_new,    # new / old
	'bin:s'=> \$bin, #可选
    'outdir=s' => \$outdir,
    'qc=s' => \$qc,                        ##yangsheng
);

#$bin ||="/NJPROJ2/HEALTH/script/Report/XKFW/V1_171011/script";   ##yangsheng
$bin ||="/PUBLIC/pipline/script/Report/XKFW/V3_180824/script";   ##sunxq
my $gxydiseasefile1="/PUBLIC/pipline/script/Report/XKFW/V3_180824/database/GXY.disease.table.3.txt";
my $ln_outdir="$outdir/03result";
mkdir "$outdir/00bin" unless -e "$outdir/00bin";
mkdir "$outdir/02individual" unless -e "$outdir/02individual";
mkdir "$outdir/03result" unless -e "$outdir/03result";

if($in_dir ne '')
{
	$right="$in_dir/report/01products/01GXY/GXY_sample_rigth_to_PDF_final.txt";
	$cusi="$in_dir/report/01products/02Cus/CS.pdf.txt";
	$DM="$in_dir/report/01products/05DM/DM.pdf.txt";
	$IDT="$in_dir/report/01products/06IDT/IDT.pdf.txt";
	$gresult="$in_dir/report/LOG/Drug.FWgxy_drug_list";
	$yresult="$in_dir/report/LOG/Drug.yaomin_list";
	$dresult="$in_dir/report/LOG/Drug.diabetes_gxy_list"
}

my (%gresultlist,%yresultlist,%dresultlist);
open GRE,$gresult;
    %gresultlist = map{chomp;my $name=basename($_);my $ID=(split/\./,$name)[0];$ID=>$_;}<GRE>;
close GRE;
open YRE,$yresult;
    %yresultlist = map{chomp;my $name=basename($_);my $ID=(split/\./,$name)[0];$ID=>$_;}<YRE>;
close YRE;
open DRE,$dresult;
    %dresultlist = map{chomp;my $name=basename($_);my $ID=(split/\./,$name)[0];$ID=>$_;}<DRE>;
close DRE;
my %yinyang_dict = ('0'=>'阴性','1'=>'阳性'); #wph add 190325

sub judge_yinyang{
    my ($ID,$right,$outDir,$type)=@_;
    open YY,"$right" || die "can't open the $right!";
    my $yinYangXing=0;
    while(<YY>){
        chomp; 
        next if ($_=~/^ID/i);
        next if($_=~/^Sample/i);
        #next if($_=~/^##/i);
        #next if (not $_  $ID);
        my @arr=split/\t/,$_;
        next unless ( $arr[0] eq $ID );
        next if ((join"",@arr[2..6]) eq '-----');
        $yinYangXing=1;
    }
    close YY;
    mkdir "$outdir/02individual/$ID" unless -e "$outdir/02individual/$ID";
    mkdir "$outdir/03result/$ID"     unless -e "$outdir/03result/$ID";
    
    open II,">$outDir/02individual/$ID/$ID.$type.yinYangXing.txt" or die "can't open the $outDir/02individual/$ID/$ID.$type.yinYangXing.txt!";
    print II "$yinYangXing\n";
    close II;
    
    return $yinYangXing;
}

my $out_file='';
my $out_file1='';
open IN,"$cfg" || die "can't open the $cfg!";
while(<IN>){
	chomp;
	next if /^\s*$/;
	next if($_=~/身份证号/);
	my ($ID,$reportflag,$name,$panel,$doctor,$salesperson)=(split/\t/,$_)[0,1,4,11,16,17]; #sxq  5.29 wph add 190116 190612 add doctor
	$reportflag=~s/阜外高血压/单基因高血压/g;
	$name=~s/\s+/-/g;
	$doctor=~s/X//g;
	$doctor=~s/\//-/g;
	$salesperson=~s/X//g;
	my $doc_sail="";
	if($doctor && $salesperson){
		$doc_sail=".$salesperson\_$doctor";
	}
	$doc_sail=""; #暂时先不出，190619
	if($_=~/糖友/){$reportflag.="+糖友";}
	if($_=~/广州市南方医院/){$reportflag.="+南方医";} #190715 wangpenghui add 
	if($_=~/吉林大学第一医院/){$reportflag.="+吉大一";} #190927 wangpenghui add
	if($_=~/山东大学齐鲁医院/ && $_=~/进院/){$reportflag.="+齐鲁";} #200409 wangpenghui add
	say '*'x40;
	say "订单${ID}所出具模块为：$reportflag";
	my($yang_gxy,$yang_cs,$yang_dm,$yang_idt,$yang_all)=(2,2,2,2,0);
	my $pedigree="";
	my @moduleall=();
	my @module_second=();
	my $module_zh="";
	my $sanger_flag=`cut -f3 $outdir/02individual/$ID/fig/$ID.family.site.txt|sort|uniq|wc -l`;#wph add 控制家系报告条件 190314
	if($panel=~/IDT/)
	{
		$yang_idt=&judge_yinyang($ID,$IDT,$outdir,'idt');
		push @moduleall,'IDT';
		open II,">$outdir/02individual/$ID/$ID.gxy.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.gxy.yinYangXing.txt!";
		print II "0\n";
		close II;

		open II,">$outdir/02individual/$ID/$ID.cs.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.cs.yinYangXing.txt!";
		print II "0\n";
		close II;
	}

	if($reportflag=~/单基因糖尿病/)
	{
		$yang_dm=&judge_yinyang($ID,$DM,$outdir,'dm');
		push @moduleall,'DM';
		open II,">$outdir/02individual/$ID/$ID.gxy.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.gxy.yinYangXing.txt!";
		print II "0\n";
		close II;

		open II,">$outdir/02individual/$ID/$ID.cs.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.cs.yinYangXing.txt!";
		print II "0\n";
		close II;
	}

	if($reportflag=~/单基因高血压三项/)
	{	
		$yang_gxy=&judge_yinyang($ID,$right,$outdir,'gxy');
		$yang_cs=&judge_yinyang($ID,$cusi,$outdir,'cs');
		push @moduleall,'All';
		push @module_second,'All';
	}
	elsif($reportflag=~/单基因高血压/ && $reportflag=~/心源性猝死/)
	{
		$yang_gxy=&judge_yinyang($ID,$right,$outdir,'gxy');
		$yang_cs=&judge_yinyang($ID,$cusi,$outdir,'cs');
		push @moduleall,'GXY';
		push @module_second,'GXY';
		push @moduleall,'CS';
		push @module_second,'CS';
	}
	elsif($reportflag=~/单基因高血压/)
	{
		$yang_gxy=&judge_yinyang($ID,$right,$outdir,'gxy');
		open II,">$outdir/02individual/$ID/$ID.cs.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.cs.yinYangXing.txt!";
		print II "0\n";
		close II;

		push @moduleall,'GXY';
		push @module_second,'GXY';
	}
	elsif($reportflag=~/心源性/ or $reportflag=~/猝死/ or $reportflag=~/心猝/)
	{
		$yang_cs=&judge_yinyang($ID,$cusi,$outdir,'cs');
		open II,">$outdir/02individual/$ID/$ID.gxy.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.gxy.yinYangXing.txt!";
                print II "0\n";
                close II;

		push @moduleall,'CS';
		push @module_second,'CS';
	}
	
	$pedigree.="-wes_file $IDT -wes_flag 1 " if($yang_idt == 1);
	$pedigree.="-dm_file $DM -dm_flag 1 " if($yang_dm  == 1);
	$pedigree.="-gxy_file $right -gxy_flag 1 " if($yang_gxy == 1);
	$pedigree.="-cs_file $cusi -cs_flag 1 " if($yang_cs == 1);
	$pedigree.="-sanger_flag 1" if($sanger_flag>2);
	
	#system "perl $bin/Cardio_01_Hyper.pl -version v2 -product gxy -site $right -disease $gxydiseasefile1 -sampleID $ID -out $outdir/02individual/$ID -part A";

	if(!$reportflag or $reportflag eq ' ' or $reportflag eq '-')
	{
		push @moduleall,'?';
		say '请检查检测项目，是否符合‘单基因高血压三项’，‘单基因高血压’，‘心源性猝死’，‘用药指导’这几项要求！';
	}
	elsif($panel=~/IDT/)
	{
		say 'IDT样本检测不包含用药模块';
		if($yang_idt==1)
		{
			$yang_all=1;
		}
	}else{
		#($reportflag=~/用药指导/) && (push @moduleall,'Druga') or ($reportflag=~/单基因高血压/ && $reportflag!~/单基因高血压三项/ && $yang_gxy==0) && (push @moduleall,'Drugb') or ($reportflag=~/心源性猝死/ && $yang_cs==0) && (push @moduleall,'Drugb') or ($reportflag=~/单基因糖尿病/ && $yang_dm==0) && (push @moduleall,'Drugb') ;
		#($reportflag=~/用药指导/) && (push @module_second,'Druga') or ($reportflag=~/单基因高血压/ && $yang_gxy==0) && (push @module_second, 'Drugb') or ($reportflag=~/心源性猝死/ && $yang_cs==0) && (push @module_second,'Drugb') or ($reportflag=~/单基因糖尿病/ && $yang_dm==0) && (push @module_second,'Drugb') if($reportflag!~/单基因高血压三项/);
			
		if($reportflag=~/用药指导/)
		{
			push @moduleall,     'Druga';
			push @module_second, 'Druga';
			if( $yang_gxy==1 || $yang_cs==1 || $yang_dm==1 ){
				$yang_all=1;
			}
		}
		else{
			if($reportflag=~/单基因高血压/ && $reportflag!~/单基因高血压三项/ && $yang_gxy==1)
			{
				$yang_all=1;
			}
			elsif($reportflag=~/心源性猝死/ && $yang_cs==1)
			{
				$yang_all=1;
			}
			elsif($reportflag=~/单基因糖尿病/ && $yang_dm==1)
			{
				$yang_all=1;
			}

			if($yang_all==0 && $reportflag!~/单基因高血压三项/ && $reportflag!~/三高用药/ &&  $reportflag!~/糖友/ && $reportflag!~/南方医/) #wph add 190321
			{ 
				push @moduleall,     'Drugb';
				push @module_second, 'Drugb';
			}
		}

	}
	my $module= join '+', @moduleall;

	my @sperl;
	if($type_new eq "new")
	{
		if($module=~m/IDT/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $IDT -sample $ID -feature $outdir/01products/feature.list -type WES -qc $in_dir/Stat_QC.xls -genelist $outdir/sample_disea_gene.txt -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.WES.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.WES.pdf $outdir/03result/$ID/$ID\_$name\_全外显子基因检测报告$doc_sail.pdf";
		}
		if($module=~m/GXY/i or $module=~m/ALL/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $right -sample $ID -feature $outdir/01products/feature.list -type GXY -qc $in_dir/Stat_QC.xls -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.GXY.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.GXY.pdf $outdir/03result/$ID/$ID\_$name\_单基因高血压基因检测报告$doc_sail.pdf";
		}
		if($module=~m/CS/i or $module=~m/ALL/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $cusi -sample $ID -feature $outdir/01products/feature.list -type CS -qc $in_dir/Stat_QC.xls -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.CS.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.CS.pdf $outdir/03result/$ID/$ID\_$name\_心源性疾病基因检测报告$doc_sail.pdf";#190506 删除阴阳性 yinyang_dict{$yang_all}
		}
		if($module=~m/DM/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $DM -sample $ID -feature $outdir/01products/feature.list -type DM -qc $in_dir/Stat_QC.xls -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.DM.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.DM.pdf $outdir/03result/$ID/$ID\_$name\_单基因糖尿病基因检测报告$doc_sail.pdf";
		}
		if($module=~m/Drug/i or $module=~m/ALL/i)
		{
			next if($panel=~/IDT/);
			push @sperl,"perl $bin/diagnosis.drug.report.pl -gresult $gresultlist{$ID} -yresult $yresultlist{$ID} -dresult $dresultlist{$ID} -sample $ID -out $outdir/02individual/$ID  > $outdir/02individual/$ID/$ID.drug.log";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.drug.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.drug.pdf $outdir/03result/$ID/$ID\_$name\_三高用药基因检测报告$doc_sail.pdf";
		}
		
		if($pedigree =~ /sanger_flag/)
		{
			push @sperl,"perl $bin/sanger_pdf.pl $pedigree -sample $ID -outdir $outdir/02individual/$ID/";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID\_JXsanger.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID\_JXsanger.pdf $outdir/03result/$ID\_$name\_一代家系验证基因检测报告$doc_sail.pdf\n";
		}
		my $script=join"\n", @sperl;
		
		$out_file.="sh $outdir/00bin/$ID.latex.sh\n";

		open SH,">$outdir/00bin/$ID.latex.sh" || die "can't open the $outdir/00bin/$ID.latex.sh!";
		print SH "$script\n";
		close SH;
		
		next;
	}

	next unless($type_new eq "old");
	my $module= join '+', @module_second;
	my $Hyper="python2 $bin/Cardio_01_GXY_v1.py  -module $module -site $right -disease $gxydiseasefile1 -sampleID $ID -out $outdir/02individual/$ID -yinyang $outdir/02individual/$ID/${ID}.gxy.yinYangXing.txt";
	my $Yaomi="perl $bin/Cardio_02_DRUG_v1.pl -module $module -gresult $gresultlist{$ID} -yresult $yresultlist{$ID} -dresult $dresultlist{$ID} -sample $ID -outdir $outdir/02individual/$ID";
	my $Cu_si="python2 $bin/Cardio_03_CS_v2.py -module $module -sample $ID  -dir $cusi -out $outdir/02individual/$ID -yinyang $outdir/02individual/$ID/${ID}.cs.yinYangXing.txt";  ##yangsheng
	my $sanger="perl $bin/Cardio_06_sanger.pl -outdir $outdir/02individual/$ID  -sample $ID -site_gxy  $right  -site_cs $cusi  ";
#	my $fendi="perl $bin/Cardio_04_fengdi.pl -sample $ID -module $module -out $outdir/02individual/$ID";  ##yangsheng

	my $count=1;
	($module=~m/GXY/i)   && ++$count && (push @sperl,"$Hyper -part $count");
	($module=~m/CS/i)    && ++$count && (push @sperl,"$Cu_si -part $count");
	($module=~m/Drug/i)  && ++$count && (push @sperl,"$Yaomi -part $count");
	if($reportflag=~/齐鲁/){--$count;} #一代模块提前 wph add 20200413
	if( $yang_gxy==1 || $yang_cs==1) {$count+=2 ; push @sperl,"$sanger -part $count -yang_gxy $yang_gxy  -yang_cs $yang_cs ";}  
#	++$count && (push @sperl,"$fendi -part $count"); ##yangsheng

	my $numbes=@sperl;
	my ($allper,$number)=($module=~m/All/i)
	   ? ( (join"\n",($Hyper." -part 2",  $Cu_si." -part 3", $Yaomi." -part 4",)), 4) ##yangsheng
	   : ( (join"\n", @sperl), $numbes+2);
	if($module=~m/All/i  && ($yang_gxy==1 || $yang_cs==1))  {$allper.="\n$sanger -part 6 -yang_gxy $yang_gxy  -yang_cs $yang_cs\n" }; 
   
	say "订单：${ID} 模块Flag：\n${module}";
	say "正在写入对应latex文件。。。。。";

	if($module=~m/All/i){$module_zh="单基因高血压三项基因检测报告";}
	elsif($module=~m/GXY/i){$module_zh="单基因高血压基因检测报告";}
	elsif($module=~m/CS/i){$module_zh="心源性疾病基因检测报告";}
	elsif($module=~m/Drug/i){$module_zh="三高用药基因检测报告";}
	else{$module_zh="$module基因检测报告";}
	my $script=<<END;
#source /ifs/TJPROJ3/HEALTH/PROJ_BK/panqi/pipeline/genetic_test/Environment.sh
#mkdir -p $outdir/02individual/$ID
#find $outdir/02individual/$ID/ -type f ! -name *.yinYangXing.txt | xargs rm -f
#rm $outdir/02individual/$ID/*

python2 $bin/Cardio_00_all_v2.py -qc $qc -cfg $cfg -feature $outdir/01products/feature.list -sample $ID -module $module -out $outdir/02individual/$ID -cs $outdir/02individual/$ID/${ID}.cs.yinYangXing.txt -gxy $outdir/02individual/$ID/${ID}.gxy.yinYangXing.txt
$allper

/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ ${ID}_Cardio_all.tex
/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ ${ID}_Cardio_all.tex
#/PUBLIC/software/public/Graphics/texlive-2014/bin/x86_64-linux/xelatex -output-directory $outdir/02individual/$ID $ID.tex

gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -sOutputFile=$outdir/02individual/$ID/${ID}_Cardio_${module}.pdf $outdir/02individual/$ID/${ID}_Cardio_all.pdf
#rm -f $ln_outdir/${ID}_Cardio_all.pdf
if [ -f $ln_outdir/$ID\_$name\_$module_zh$doc_sail.pdf ];then
rm $ln_outdir/$ID\_$name\_$module_zh$doc_sail.pdf
fi
ln -s $outdir/02individual/$ID/${ID}_Cardio_${module}.pdf $ln_outdir/$ID\_$name\_$module_zh$doc_sail.pdf
END
	if(not $module){$script="";}
	$module= join '+', @moduleall;
	if($module=~m/DM/i)
	{
		$script.="perl $bin/diagnosis.report.pl -site $DM -sample $ID -feature $outdir/01products/feature.list -type DM -qc $in_dir/Stat_QC.xls -hospital $bin/product.config.txt -out $outdir/02individual/$ID\n";
		$script.="/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.DM.tex\n";
		$script.="ln -s $outdir/02individual/$ID/$ID.DM.pdf $outdir/03result/$ID\_$name\_单基因糖尿病基因检测报告$doc_sail.pdf\n";
	}

	if($reportflag=~/齐鲁/){
		$script=~s/Cardio_00_all_v2.py/Cardio_00_all_v2.qilu.py/g;

	}

	if($reportflag=~/南方/){
		if($module=~m/GXY/i){
			$script.="python3 $bin/Cardio_PMC_WES_dox.py -module GXY -sample $ID -dir $right  -out $outdir/02individual/$ID -yinyang $outdir/02individual/$ID/$ID.gxy.yinYangXing.txt  -qc $in_dir/Stat_QC.xls -cfg $outdir/01products/cfg.txt -feature $outdir/01products/feature.list -genelist /PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/gene_deseases/191016_v7_gene_desease\n";
			$script.="ln -s $outdir/02individual/$ID/$ID.单基因高血压-血钾异常Panel报告.docx $outdir/03result/南方医院-诺禾心康-$ID\_$name\_单基因高血压血钾异常基因检测报告.docx\n";	
		}
		if($module=~m/CS/i){
			$script.="python3 $bin/Cardio_PMC_WES_dox.py -module CS -sample $ID -dir $cusi -out $outdir/02individual/$ID -yinyang $outdir/02individual/$ID/$ID.cs.yinYangXing.txt -qc $in_dir/Stat_QC.xls -cfg $outdir/01products/cfg.txt -feature $outdir/01products/feature.list -genelist /PUBLIC/pipline/database/siteFilterDB/product/CuS/cusi_xinkang.geneV4\n";
			$script.="ln -s $outdir/02individual/$ID/$ID\.遗传性心血管基因检测Panel报告.docx $outdir/03result/南方医院-诺禾心康-$ID\_$name\_遗传性心血管基因检测报告.docx\n";
		}
	}
	if($reportflag=~/吉大一/ && $module=~m/Drug/i){ #wph add 190927
		$script.="python3 $bin/JiDaYiDrug.py -dir $gresultlist{$ID} -out $outdir/02individual/$ID -sample $ID -cfg $outdir/01products/cfg.txt\n";
		$script.="ln -s $outdir/02individual/$ID/$ID-吉大一用药报告.docx $outdir/03result/$ID-吉大一用药报告.docx\n";
	}

	if($module=~m/IDT/i)
	{   # wangpenghui add 20180531  #190314 wph 删除其他模块报告，只留全外
		$script="perl $bin/diagnosis.report.pl -site $IDT -sample $ID -feature $outdir/01products/feature.list -type WES -qc $in_dir/Stat_QC.xls -genelist $outdir/sample_disea_gene.txt -hospital $bin/product.config.txt -out $outdir/02individual/$ID\n";
		$script.="/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.WES.tex\n";
		$script.="ln -s $outdir/02individual/$ID/$ID.WES.pdf $outdir/03result/$ID\_$name\_全外显子基因检测报告$doc_sail.pdf\n";
		if($reportflag=~/南方/){
			$script.="python3 $bin/Cardio_PMC_WES_dox.py -module IDT -sample $ID -dir $IDT -out $outdir/02individual/$ID -yinyang $outdir/02individual/$ID.idt.yinYangXing.txt -qc $in_dir/Stat_QC.xls -cfg $outdir/01products/cfg.txt -feature $outdir/01products/feature.list -genelist $outdir/sample_disea_gene.txt\n";
			$script.="ln -s $outdir/02individual/$ID/$ID\.遗传病性疾病全外显子基因检测报告.docx $outdir/03result/南方医院-诺禾心康-$ID\_$name\_遗传病性疾病全外显子基因检测报告.docx\n";
		}
	}
	
	if($pedigree =~/sanger_flag/)
	{
		$script.="perl $bin/sanger_pdf.pl $pedigree -sample $ID -outdir $outdir/02individual/$ID\n";
		$script.="/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID\_JXsanger.tex\n";
		$script.="ln -s $outdir/02individual/$ID/$ID\_JXsanger.pdf $outdir/03result/$ID\_$name\_家系一代验证基因检测报告$doc_sail.pdf\n";
	}

	##yangsheng '-qc'
	open SHH,">$outdir/00bin/$ID.latex.1.sh" || die "can't open the $outdir/00bin/$ID.latex.1.sh!";
	print SHH "$script";
	close SHH;
	
	$out_file1.="sh $outdir/00bin/$ID.latex.1.sh\n";
}
close IN;

if($out_file1 ne '')
{
	open FILE, "> $outdir/00bin/latex.1.pipline.sh" || die "can't open the $outdir/00bin/latex.1.pipline.sh!";
	print FILE $out_file1;
	close FILE;
}

if($out_file ne '')
{
	open OUT,  "> $outdir/00bin/latex.pipline.sh"   || die "can't open the $outdir/00bin/latex.pipline.sh!";
	print OUT $out_file;
	close OUT;
}

say '*'x40;
say "latex文件写入完成！！";
