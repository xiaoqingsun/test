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
	$right="$outdir/01products/01GXY/GXY_sample_rigth_to_PDF_final.txt";
	$cusi="$outdir/01products/02Cus/CS.pdf.txt";
	$DM="$outdir/01products/05DM/DM.pdf.txt";
	$IDT="$outdir/01products/06IDT/IDT.pdf.txt";
	$gresult="$outdir/LOG/Drug.FWgxy_drug_list";
	$yresult="$outdir/LOG/Drug.yaomin_list";
	$dresult="$outdir/LOG/Drug.diabetes_gxy_list"
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
	my ($ID,$reportflag,$panel)=(split/\t/,$_)[0,1,11]; #sxq  5.29
	$reportflag=~s/阜外高血压/单基因高血压/g;
	say '*'x40;
	say "订单${ID}所出具模块为：$reportflag";
	my($yang_gxy,$yang_cs,$yang_dm,$yang_idt)=(2,2,2,2);
	my $pedigree="";
	my @moduleall=();
	my @module_second=();
	
	if($panel=~/IDT/)
	{
		$yang_idt=&judge_yinyang($ID,$IDT,$outdir,'idt');
		push @moduleall,'IDT';
                push @module_second,'WES';
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
                push @module_second,'DM';
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
		push @module_second,'GXY';
                push @module_second,'CS';
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
	
	#system "perl $bin/Cardio_01_Hyper.pl -version v2 -product gxy -site $right -disease $gxydiseasefile1 -sampleID $ID -out $outdir/02individual/$ID -part A";

	if(!$reportflag or $reportflag eq ' ' or $reportflag eq '-')
	{
		push @moduleall,'?';
		say '请检查检测项目，是否符合‘单基因高血压三项’，‘单基因高血压’，‘心源性猝死’，‘用药指导’这几项要求！';
	}
	elsif($panel=~/IDT/)
	{
		say 'IDT样本检测不包含用药模块';
	}else{
		($reportflag=~/用药指导/) && (push @moduleall,'Druga') or ($reportflag=~/单基因高血压/ && $reportflag!~/单基因高血压三项/ && $yang_gxy==0) && (push @moduleall,'Drugb') or ($reportflag=~/心源性猝死/ && $yang_cs==0) && (push @moduleall,'Drugb') or ($reportflag=~/单基因糖尿病/ && $yang_dm==0) && (push @moduleall,'Drugb') ;
                if($pedigree eq '' && $reportflag!~/用药指导/){
			push @module_second,'Drugb';
		}elsif($reportflag=~/用药指导/){
			push @module_second,'Druga';
		}
	}

	my @sperl;
	if($type_new eq "new")
	{
		my $module= join '+', @moduleall;
		if($module=~m/IDT/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $IDT -sample $ID -type WES -qc $in_dir/Stat_QC.xls -genelist $outdir/sample_disea_gene.txt -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.WES.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.WES.pdf $outdir/03result/$ID/$ID.WES.pdf";
		}
		if($module=~m/GXY/i or $module=~m/ALL/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $right -sample $ID -type GXY -qc $in_dir/Stat_QC.xls -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.GXY.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.GXY.pdf $outdir/03result/$ID/$ID.GXY.pdf";
		}
		if($module=~m/CS/i or $module=~m/ALL/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $cusi -sample $ID -type CS -qc $in_dir/Stat_QC.xls -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.CS.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.CS.pdf $outdir/03result/$ID/$ID.CS.pdf";
		}
		if($module=~m/DM/i)
		{
			push @sperl,"perl $bin/diagnosis.report.pl -site $DM -sample $ID -type DM -qc $in_dir/Stat_QC.xls -hospital $bin/product.config.txt -out $outdir/02individual/$ID";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.DM.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.DM.pdf $outdir/03result/$ID/$ID.DM.pdf";
		}
		if($module=~m/Drug/i or $module=~m/ALL/i)
		{
			next if($panel=~/IDT/);
			push @sperl,"perl $bin/diagnosis.drug.report.pl -gresult $gresultlist{$ID} -yresult $yresultlist{$ID} -dresult $dresultlist{$ID} -sample $ID -out $outdir/02individual/$ID  > $outdir/02individual/$ID/$ID.drug.log";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID.drug.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID.drug.pdf $outdir/03result/$ID/$ID.drug.pdf";
		}
		
		if($pedigree ne "")
		{
			push @sperl,"perl $bin/sanger_pdf.pl $pedigree -sample $ID -outdir $outdir/02individual/$ID/";
			push @sperl,"/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID\_JXsanger.tex";
			push @sperl,"ln -s $outdir/02individual/$ID/$ID\_JXsanger.pdf $outdir/03result/$ID\_JXsanger.pdf\n";
		}
		my $script=join"\n", @sperl;
		
		$out_file.="sh $outdir/00bin/$ID.latex.sh\n";

		open SH,">$outdir/00bin/$ID.latex.sh" || die "can't open the $outdir/00bin/$ID.latex.sh!";
		print SH "$script\n";
		close SH;
		
		next;
	}

	next unless($type_new eq "old");
	if(!-e "$outdir/02individual/$ID/$ID.dm.yinYangXing.txt"){
		open II,">$outdir/02individual/$ID/$ID.dm.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.dm.yinYangXing.txt!";
	        print II "0\n";
        	close II;
	}
	if(!-e "$outdir/02individual/$ID/$ID.wes.yinYangXing.txt"){
		open II,">$outdir/02individual/$ID/$ID.wes.yinYangXing.txt" or die "can't open the $outdir/02individual/$ID/$ID.wes.yinYangXing.txt!";
	        print II "0\n";
        	close II;
	}
	my $module= join '+', @module_second;
	my $Hyper="python2 $bin/Cardio_01_GXY_v1.py  -module $module -site $right -disease $gxydiseasefile1 -sampleID $ID -out $outdir/02individual/$ID -yinyang $outdir/02individual/$ID/${ID}.gxy.yinYangXing.txt";
	my $Yaomi="perl $bin/Cardio_02_DRUG_v1.pl -module $module -gresult $gresultlist{$ID} -yresult $yresultlist{$ID} -dresult $dresultlist{$ID} -sample $ID -outdir $outdir/02individual/$ID";
	my $Cu_si="python2 $bin/Cardio_03_CS_v2.py -module $module -sample $ID  -dir $cusi -out $outdir/02individual/$ID -yinyang $outdir/02individual/$ID/${ID}.cs.yinYangXing.txt";  ##yangsheng
        my $diabetes="python2 $bin/Cardio_04_DM_v1.py  -module  $module  -dir $DM   -sample $ID  -out  $outdir/02individual/$ID  -yinyang  $outdir/02individual/$ID/${ID}.dm.yinYangXing.txt"; 
        my $wesshell="python2 $bin/Cardio_05_WES_v1.py  -module  $module  -dir $IDT   -sample $ID  -out  $outdir/02individual/$ID  -yinyang  $outdir/02individual/$ID/${ID}.idt.yinYangXing.txt";
	my $sanger="perl $bin/Cardio_06_sanger.pl -outdir $outdir/02individual/$ID  -sample $ID -site_gxy  $right  -site_cs $cusi    -site_dm  $DM  -site_wes $IDT  ";
#	my $fendi="perl $bin/Cardio_04_fengdi.pl -sample $ID -module $module -out $outdir/02individual/$ID";  ##yangsheng

	my $count=1;
        if($module=~m/All/i){
		++$count && (push @sperl,"$Hyper -part $count");
                ++$count && (push @sperl,"$Cu_si -part $count");
                ($module=~m/DM/i)   && ++$count && (push @sperl,"$diabetes -part $count");
        	($module=~m/IDT/i)  && ++$count && (push @sperl,"$wesshell -part $count") && (print "注意IDT与三项的用药模块");	
                ($module=~m/Drug/i)  && ++$count && (push @sperl,"$Yaomi -part $count");
        }else{
		($module=~m/GXY/i)   && ++$count && (push @sperl,"$Hyper -part $count");
	        ($module=~m/CS/i)    && ++$count && (push @sperl,"$Cu_si -part $count");
                ($module=~m/DM/i)   && ++$count && (push @sperl,"$diabetes -part $count");
                ($module=~m/IDT/i)    && ++$count && (push @sperl,"$wesshell -part $count");
                ($module=~m/Drug/i)  && ++$count && (push @sperl,"$Yaomi -part $count");
	}
	if( $pedigree ne '') {$count+=2 ; push @sperl,"$sanger -part $count -yang_gxy $yang_gxy  -yang_cs $yang_cs -yang_dm $yang_dm  -yang_wes $yang_idt ";}  
  	my $allper=join"\n", @sperl; 
	say "订单：${ID} 模块Flag：\n${module}";
	say "正在写入对应latex文件。。。。。";

	my $script=<<END;
#source /ifs/TJPROJ3/HEALTH/PROJ_BK/panqi/pipeline/genetic_test/Environment.sh
mkdir -p $outdir/02individual/$ID
find $outdir/02individual/$ID/ -type f ! -name *.yinYangXing.txt | xargs rm -f
#rm $outdir/02individual/$ID/*

python2 $bin/Cardio_00_all_v2.py -qc $qc -cfg $cfg -sample $ID -module $module -out $outdir/02individual/$ID -cs $outdir/02individual/$ID/${ID}.cs.yinYangXing.txt -gxy $outdir/02individual/$ID/${ID}.gxy.yinYangXing.txt  -dm $outdir/02individual/$ID/${ID}.dm.yinYangXing.txt -wes $outdir/02individual/$ID/${ID}.wes.yinYangXing.txt
$allper

/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ ${ID}_Cardio_all.tex
/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ ${ID}_Cardio_all.tex
#/PUBLIC/software/public/Graphics/texlive-2014/bin/x86_64-linux/xelatex -output-directory $outdir/02individual/$ID $ID.tex

gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -sOutputFile=$outdir/02individual/$ID/${ID}_Cardio_${module}.pdf $outdir/02individual/$ID/${ID}_Cardio_all.pdf
#rm -f $ln_outdir/${ID}_Cardio_all.pdf
ln -s $outdir/02individual/$ID/${ID}_Cardio_${module}.pdf $ln_outdir/
END
	
	
	if($pedigree ne "")
	{
		$script.="perl $bin/sanger_pdf.pl $pedigree -sample $ID -outdir $outdir/02individual/$ID\n";
		$script.="/PUBLIC/pipline/software/texlive-2014/bin/x86_64-linux/xelatex -interaction=nonstopmode -output-directory $outdir/02individual/$ID/ $ID\_JXsanger.tex\n";
		$script.="ln -s $outdir/02individual/$ID/$ID\_JXsanger.pdf $outdir/03result/$ID\_JXsanger.pdf\n";
	}

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
