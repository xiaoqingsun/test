#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use autodie;
use File::Basename;
use File::Path qw(make_path);
my ($indir,$outdir,$bin,$types,$list,$switch);
GetOptions(
	'indir=s'  => \$indir,
	'outdir=s' => \$outdir,
	'type=s' => \$types,
	'bin=s'	=>	\$bin,
	'list=s' =>	\$list,
	'step=i' => \$switch,
);
#$bin ||="/ifs/TJPROJ3/HEALTH/zhongying/01script/06pipline/04FYGXY/report/script/";
$bin ||="/PUBLIC/pipline/script/Report/XKFW/V3_180824/script/"; #zhouwei change 20170711 ##yangsheng
$list ||="all";
my $typesall="";
$switch ||=2;
chomp($outdir=`pwd`) unless ($outdir);
$outdir =~ s/\/$//;
chomp(my $cwd =`pwd`);
$outdir="$cwd/$outdir" unless($outdir=~/^\//);

my $usage=<<END;
perl $0 
	-indir	输入目录，多个用 ',' 分隔
	-outdir 报告输出目录
	-type the type of products, GXY/CS/Drug 写任意产品组合，用，分割开
	-bin the pathway of the script 
	-list  all或者不写该参数 表示 输出该路径下所有的样本；
		   输入list文件， 表示 输出list文件内指定的样本；
	-switch  默认为2，即输出PDF+微信端，
			 如果参数为1，即只输出需要解读组check的点；
END

#=================judge argument ===================
die "$usage\n\t\t-indir must be give " unless $indir;
die "$usage\n\t\t-outdir must be give " unless $types;

#================prepare pathway =============================
if($switch==1){
	make_path("$outdir/LOG","$outdir/00bin","$outdir/01products","$outdir/check","$outdir/check/HgmdClivar_site");
}else{
	make_path("$outdir/LOG","$outdir/00bin","$outdir/01products","$outdir/check");
};

# zhouwei change 20170712

my $step01cs="$bin/step1_get_truesite_and_CS_up.pl";
my $step01cs_2="$bin/step2_get_truesite_and_CSDM_up.pl";

my $step01wes="$bin/step1_get_truesite_and_WES_up.pl";
my $step02wes="$bin/step2_get_truesite_and_WES_up.pl";
my $step01filter_lowdp="$bin/filter_low_dp.sh";#wang penghui change 20181227
my $step01addgene="$bin/add_gene_sample.py"; #wang penghui change 20180531
my $step01addtype="$bin/add_sample_panel.py"; #wang penghui change 20181112

my $step01dm="$bin/step1_get_truesite_and_DM_up.pl";
my $step02dm="$bin/step2_get_truesite_and_CSDM_up.pl";

my $fw_xk_gxy="$bin/Generate_XK_FW_Hypertation_PDFV2.pl";
my $refer="$bin/reference.pl";
my $txt2xlsx="$bin/txt2xlsx.py";

my %hash=(
	"CS"=>"04.Annotation",
	"GXY"=>"04.Annotation",
	"Drug"=>"05.Drug",
	"DM"=>"04.Annotation",
    "IDT"=>"04.Annotation",
);
my %mkpath=(
	"GXY"=>"01GXY",
	"CS"=>"02Cus",
	"Drug"=>"03Drug",
    "DM"=>"05DM",
    "IDT"=>"06IDT",
);
########################## main #####################################################
my @files;
my %samples;
my %samplesDM;
my @final_2_txt;
open BIN," > $outdir/00bin/pipline.sh" || die "can't open the file!";
for ( (split /,/,$indir)){
	print glob "$_"."\n";
	print "$indir\n";
	push @files , glob "$_"; #修改此处来控制文件输入类型
}
foreach my $n(@files){
	my @files=`find $n/cr/ -name  "*.ALL.final.info.opt" `;
	push @final_2_txt,@files;
}

&sample_list(@files);
&cfg(@files);
if($types =~/IDT/){
	&script("CS") if($typesall=~/CS/);
	&WES_script("IDT");
	&GXY_script("GXY") if($typesall=~/GXY/);
	&DM_script("DM") if($typesall=~/DM/);
}else{
	&script("CS");
	&GXY_script("GXY");
	&DM_script("DM");
}
&file_list("Drug",@files);
&file_list("DM",@files);
if($switch ==1){
`iconv -f UTF-8 -t gbk -c $outdir/01products/cfg.txt > $outdir/check/cfg.xls`;
`ln -s $outdir/../cr/N*/*/04.Annotation/*.ALL.final.info.opt.HgmdClivar.xls $outdir/check/HgmdClivar_site/`;
}#wph add 20181115
if($switch !=1){
	open PDF, "> $outdir/00bin/producePDF.sh" || die "can't open the $outdir/00bin/producePDF.sh !";
	print PDF <<END;
perl $bin/Cardio_00_all_feature.pl $indir/report/01products/cfg.txt v2 1> $indir/report/01products/feature.list 2> $indir/report/01products/feature.list.err
perl $bin/Cardio_05_All_v2.pl -cfg $outdir/01products/cfg.txt -qc $indir/Stat_QC.xls -in $indir -type new -module $types -outdir $outdir
perl $bin/Cardio_05_All_v2.pl -cfg $outdir/01products/cfg.txt -qc $indir/Stat_QC.xls -in $indir -type old -module $types -outdir $outdir
END

#perl $bin/Cardio_05_All_v1.pl -cfg $indir/report/01products/cfg.txt -qc $indir/Stat_QC.xls -right $indir/report/01products/01GXY/GXY_sample_rigth_to_PDF_final.txt  -gresult $outdir/LOG/Drug.FWgxy_drug_list -yresult $outdir/LOG/Drug.yaomin_list -dresult $outdir/LOG/Drug.diabetes_gxy_list -cusi $indir/report/01products/02Cus/CS.pdf.txt -module $types -outdir $outdir

##yangsheng '-qc'
	close PDF;
	print BIN "sh $outdir/00bin/producePDF.sh\nsh $outdir/00bin/latex.pipline.sh\n"; #zhouwei change 20170410
}

print BIN "#perl /PUBLIC/pipline/script/Report/XKFW/V3_180824/script/substr_fa_path.pl  $outdir/\n";#wph dell 190319
print BIN "sh $step01filter_lowdp $indir\n" ;#wph add 190904
close BIN;
#system "nohup sh $outdir/00bin/pipline.sh &";

##################################### sub ######################################
sub script{
	my $product=shift;
	mkdir "$outdir/01products/$mkpath{$product}" unless -e "$outdir/01products/$mkpath{$product}";
	&file_list($product,@files); # 防止参数传递错误，需要数组在后面
	open YC,">$outdir/00bin/$mkpath{$product}.site.sh";
	my $yc_1=<<END;
source /PUBLIC/pipline/Environment.sh
perl $step01cs -path $outdir/LOG/$product.file_list  -type $product -list  $outdir/LOG/sample_list -out $outdir/01products/$mkpath{$product}
python $step01addtype -t $product -p $outdir/01products/$mkpath{$product}/
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/$product.right.txt -t $outdir/check/$product.right.xlsx
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/$product.err.txt -t $outdir/check/$product.err.xlsx
#iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.err.txt > $outdir/check/$product.err.xls
#iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.right.txt >$outdir/check/$product.right.xls
END
    my $yc_2=<<END;
source /PUBLIC/pipline/Environment.sh
perl $step01cs_2 -type $product -out $outdir/01products/$mkpath{$product}
if [ `cut -f83 $outdir/01products/$mkpath{$product}/$product.pdf.txt |grep -v '-'|grep -c '[0-9]'` -gt 0 ];then
perl $refer      -type $product -out $outdir/01products/$mkpath{$product}  -reference /PUBLIC/work/sunxiaoqing/pubmed/pubmed_20180928_reference 
fi
iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.pdf.txt >$outdir/check/$product.pdf.xls
END
	if($switch==1){
		print YC "$yc_1";
	}else{
		print YC "$yc_2";
	}
	close YC;
	print BIN "sh $outdir/00bin/$mkpath{$product}.site.sh\n";
}

sub WES_script{
	my $product=shift;
        mkdir "$outdir/01products/$mkpath{$product}" unless -e "$outdir/01products/$mkpath{$product}";
        &file_list($product,@files);
        open YC,">$outdir/00bin/$mkpath{$product}.site.sh";	
        my $wes_1=<<END;
source /PUBLIC/pipline/Environment.sh
python $step01addgene $outdir/01products/cfg.txt $outdir/sample_disea_gene.txt 
perl $step01wes -path $outdir/LOG/$product.file_list  -type $product -out $outdir/01products/$mkpath{$product} -samplelist  $outdir/LOG/sample_list  -genelist $outdir/sample_disea_gene.txt 
python $step01addtype -t $product -p $outdir/01products/$mkpath{$product}/
#wang penghui change 20180531
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/$product.err.txt -t $outdir/check/$product.err.xlsx
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/$product.right.txt -t $outdir/check/$product.right.xlsx
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/$product.right.2.txt -t $outdir/check/$product.right.all.xlsx
#iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.err.txt > $outdir/check/$product.err.xls
#iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.right.txt >$outdir/check/$product.right.xls
#iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.right.2.txt >$outdir/check/$product.right.all.xls
END
    #wph add 180821
    #wph changed product 2 typesall 190328
	my $wes_2=<<END;
source /PUBLIC/pipline/Environment.sh
perl $step02wes -type $typesall$product -out $outdir/01products/$mkpath{$product}
if [ `cut -f83 $outdir/01products/$mkpath{$product}/$product.pdf.txt |grep -v '-'|grep -c '[0-9]'` -gt 0 ];then
perl $refer      -type $product -out $outdir/01products/$mkpath{$product}  -reference /PUBLIC/work/sunxiaoqing/pubmed/pubmed_20180928_reference
fi
iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.pdf.txt >$outdir/check/$product.pdf.xls
END
	if($switch==1){
        print YC "$wes_1";
    }else{
		print YC "$wes_2";
	}
    close YC;
    print BIN "sh $outdir/00bin/$mkpath{$product}.site.sh\n";
}

sub GXY_script{
	my $product=shift;	
	mkdir "$outdir/01products/$mkpath{$product}" unless -e "$outdir/01products/$mkpath{$product}";
	&file_list($product,@files);
	open GXY,">$outdir/00bin/$mkpath{$product}.site.sh";
	my $gxy_1=<<END;
source /PUBLIC/pipline/Environment.sh
perl $fw_xk_gxy -final_2_list $outdir/LOG/$product.file_list -out $outdir/01products/$mkpath{$product} -step $switch
python $step01addtype -t $product -p $outdir/01products/$mkpath{$product}/
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/GXY_sample_rigth_info.txt -t $outdir/check/GXY_sample_rigth_info.xlsx
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/GXY_sample_wrong_info.txt -t $outdir/check/GXY_sample_wrong_info.xlsx
#iconv -f UTF-8 -t gbk -c $outdir/01products/$mkpath{$product}/GXY_sample_rigth_info.txt > $outdir/check/GXY_sample_rigth_info.xls
#iconv -f UTF-8 -t gbk -c $outdir/01products/$mkpath{$product}/GXY_sample_wrong_info.txt > $outdir/check/GXY_sample_wrong_info.xls
END
	my $gxy_2=<<END;
source /PUBLIC/pipline/Environment.sh
perl $fw_xk_gxy -final_2_list $outdir/LOG/$product.file_list -out $outdir/01products/$mkpath{$product} -step $switch
if [ `cut -f83 $outdir/01products/$mkpath{$product}/GXY_sample_rigth_to_PDF_final.txt |grep -v '-'|grep -c '[0-9]'` -gt 0 ];then
perl $refer      -type $product -out $outdir/01products/$mkpath{$product}  -reference /PUBLIC/work/sunxiaoqing/pubmed/pubmed_20180928_reference
fi
#iconv -f UTF-8 -t gbk -c $outdir/01products/$mkpath{$product}/GXY_sample_rigth_core_site.txt > $outdir/check/GXY_sample_rigth_core_site.xls
iconv -f UTF-8 -t gbk -c $outdir/01products/$mkpath{$product}/GXY_sample_rigth_to_PDF_final.txt > $outdir/check/GXY_sample_rigth_to_PDF_final.xls
END
##zhouwei changge 20170607
	if($switch==1){
		print GXY "$gxy_1";
	}else{
		print GXY "$gxy_2";
	}
	close GXY;
	print BIN "sh $outdir/00bin/$mkpath{$product}.site.sh\n";
}

sub DM_script{
	my $product=shift;
	if($switch==1 && -e "$outdir/01products/cfgDM.txt"){
		mkdir "$outdir/01products/$mkpath{$product}" unless -e "$outdir/01products/$mkpath{$product}";
		open DM,">$outdir/00bin/$mkpath{$product}.site.sh";
       		my $dm_1=<<END;
source /PUBLIC/pipline/Environment.sh
perl $step01dm -path $outdir/LOG/$product.file_list  -type $product -list  $outdir/LOG/sample_list_DM -out $outdir/01products/$mkpath{$product}
python $step01addtype -t $product -p $outdir/01products/$mkpath{$product}/
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/$product.err.txt -t $outdir/check/$product.err.xlsx
python2 $txt2xlsx $outdir/01products/$mkpath{$product}/$product.right.txt -t $outdir/check/$product.right.xlsx
#iconv -f UTF-8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.err.txt > $outdir/check/$product.err.xls
#iconv -f UTF-8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.right.txt >$outdir/check/$product.right.xls
END
        	print DM "$dm_1";
        	close  DM;
	}elsif($switch==2  && -e "$outdir/01products/cfgDM.txt"){
        	open DM,">$outdir/00bin/$mkpath{$product}.site.sh";  
                my $dm_2=<<END;
source /PUBLIC/pipline/Environment.sh
perl $step02dm -type $product -out $outdir/01products/$mkpath{$product}
if [ `cut -f83 $outdir/01products/$mkpath{$product}/$product.pdf.txt |grep -v '-'|grep -c '[0-9]'` -gt 0 ];then
perl $refer      -type $product -out $outdir/01products/$mkpath{$product}  -reference /PUBLIC/work/sunxiaoqing/pubmed/pubmed_20180928_reference
fi
iconv -f utf8 -t gbk -c $outdir/01products/$mkpath{$product}/$product.pdf.txt >$outdir/check/$product.pdf.xls
END
		print DM "$dm_2";
	}
        print BIN "sh $outdir/00bin/$mkpath{$product}.site.sh\n"; 
}

########################## deal CFG ###################################
sub cfg{
	my @files=@_; 
	my $num=0;
        my $numDM=0;
	my @jnlp;
	open CFG,"> $outdir/01products/cfg.txt" || die "can't open the $outdir!";
	foreach my $n(@files){ 
		if (-s "$n/result/cfg_add_sex.xls"){
			push @jnlp , glob "$n/result/*jnlp";
			open IN,"$n/result/cfg_add_sex.xls";
			while(<IN>){
				if($num==0 && $_=~/身份证号/){
					print CFG "$_";
				}else{
					next if($_=~/身份证号/);
					print CFG "$_";
				}
                $numDM++ if($_=~/单基因糖尿病/);
				if($_=~/单基因糖尿病/){$typesall.="DM+";}
				if($_=~/单基因高血压/){$typesall.="GXY+";}
				if($_=~/心源性猝死/){$typesall.="CS+";}
			}
			close IN;
			$num++;
		}else{
			print "please check the $n, whether exists the result/cfg_add_sex.xls; if not, you can run the stat.QC.sh !";
		}
	}
	close OUT;
	close CFG;
	for my $i(@jnlp){
		my $ID=(split/\./,basename($i))[0];
		if(exists $samples{$ID}){
			`cp $i $outdir/check/$ID.jnlp`;
		}else{
			next;
		}
	}
        if($numDM>0){
                `grep '单基因糖尿病' $outdir/01products/cfg.txt > $outdir/01products/cfgDM.txt`;
		`cut -f1 $outdir/01products/cfgDM.txt > $outdir/LOG/sample_list_DM`;
                 open CONF,"$outdir/LOG/sample_list_DM" || die "can't open the $outdir/LOG/sample_list_DM!";
			%samplesDM=map{chomp; (split/\t/,$_)[0] => 1;}<CONF>;
    		  close CONF;
	}
}

################################### get the files for stat ########################
sub file_list{
	my ($type,@files)=@_;
	if($type eq "Drug"){
		open OUT,"> $outdir/LOG/$type.yaomin_list" || die "can't open the file!";
		open PDF,"> $outdir/LOG/$type.FWgxy_drug_list" || die "can't open the file!";
		open OUT2,"> $outdir/LOG/$type.diabetes_gxy_list" || die "can't open the file!";
		foreach my $n(@files){
			print "##$n\n";
			my @yaomin_files=`find $n/cr/ -name "*yaomin_gxy.txt"`;
			my @Fw_files=`find $n/cr/ -name "*FWgxy_drug.pdf.txt"`;
			my @dia_files=`find $n/cr/ -name "*diabetes_gxy.txt"`;
			my @finals=&file_filter($type,\@yaomin_files,\%samples);
			my @pdfs=&file_filter($type,\@Fw_files,\%samples);
			my @dias=&file_filter($type,\@dia_files,\%samples);
			print OUT join("\n",@finals)."\n";
			print PDF join("\n",@pdfs)."\n";
			print OUT2 join("\n",@dias)."\n";
		}
		close OUT;
		close PDF;
		close OUT2;	
	}elsif($type =~/GXY/){
		open OUT,"> $outdir/LOG/$type.file_list" || die "can't open the file!";
		my @finals=&file_filter($type,\@final_2_txt,\%samples); 
		print OUT join("\n",@finals);
		close OUT;
	}elsif($type =~/DM/){
		open OUT,"> $outdir/LOG/$type.file_list" || die "can't open the file!";
		foreach my $key(keys %samplesDM){print "$key\n";}
                my @finals=&file_filter($type,\@final_2_txt,\%samplesDM);
                print OUT join("\n",@finals);
                close OUT;
	}else{
		open OUT,"> $outdir/LOG/$type.file_list" || die "can't open the file!";
		my @finals=&file_filter($type,\@final_2_txt,\%samples);
		print OUT join("\n",@finals)."\n";
		close OUT;
	}
}

sub file_filter{
	my ($type,$file,$samp)=@_;
	my @right;
	for my $k(@$file){
		chomp $k;
		if($k=~/$hash{$type}/i){
			my $ID=(split/\./,basename($k))[0]; 
			if(exists $$samp{$ID}){
				push @right,$k;  
			}else{
				next;
			}
		}else{
			next;
		}
	}
	return @right;
}

############################ get the sample for weixinduan#########################
sub sample_list{
	my @files=@_;
	my @temp;
	foreach my $n(@files){
		push @temp , glob "$n/cr/*";
	};
	my @arr=grep !/ReportForge/,@temp;
	open OUT,"> $outdir/LOG/sample_list" || die "can't open the file!";
	if($list ne "all"){
		my $n;
		open IN,"$list" || die "can't open the file!";
		while(<IN>){
			chomp;
			for my $i(@arr){
				if($i=~/$_/){
					$n++;
					last;
				}else{
					next;
				}
			}
			if($n!=0){
				$samples{$_}=1;
				print OUT "$_\n";
			}else{
				print "$_ is not enclosed the path, please to check!\n";
			}
		}
		close IN;
	}else{
		for my $i(@arr){
			my $ID=basename($i);
			$samples{$ID}=1;
			print OUT "$ID\n";
		}
	}
}
