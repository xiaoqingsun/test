use strict;
use warnings;
#use utf8;
#use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util qw/max min sum maxstr minstr shuffle/;
use feature qw/say/;

my ($path,$out,$types,$cfg,$genelist,$samplelist);
GetOptions(
	'path=s'  => \$path,
	'out=s' => \$out,
	'type=s'	=> \$types,
	'cfg=s' => \$cfg,
        'genelist=s' => \$genelist,
	'samplelist=s' => \$samplelist,
);

my $usage=<<END;

perl $0 
	-path	需要检测样本的路径list文件，例如/ifs/TJPROJ3/HEALTH/Project/2C/160926_ST-E00192_0305_BH35GVALXX/cr/NHE0000000649/NHE0000000649/04YC/NHE0000000649.final.2.txt
	-type	产品类型,YC | YC133(遗传病)、GXY、CS、CR 、NNY
	-out	输出路径,输出相应产品的true的点，且YC，GXY，CS会与最新位点数据库匹配,返回与数据库匹配的结果
END

my %disease;
my %database;
my %oldbase;
my %fubiaobasee;
my %database_2;
my %sub_sampleline;
my %case_ditacion_ratio;
my %RTsampleIDs;
my %dupsites;
my %product=(
	"CS"=>"CS",
	"GXY"=>"GXY",
	"IDT"=>"IDT",
);
my $type=$product{$types};
my $disease_file="/PUBLIC/pipline/database/data/gene/gene_add_desase.xls";  #zhouwei change 20170712
my $control_detetion="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/control_detection_ratio/170122_v0_1256wes_detection_ratio_out"; #zhouwei change 20170712 
my $case_detection ||="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/gxy_detection_ratio/";    #case集样本位点检出率
my $cusigene=$genelist;
my %gene=&genelist($cusigene);
my %detetion=&Deteion($control_detetion);
&sub_case_ditacion_ratio;
my %samplelist=&readsamplelist($samplelist);

################################################# main #######################################################################
open OUT, ">$out/$types.err.txt" || die "can't open the $out/$types.err.txt!";
open RT, "> $out/$types.right.txt" || die "can't open the $out/$types.right.txt!";
open RT2, "> $out/$types.right.2.txt" || die "can't open the $out/$types.right.2.txt!";#add right.2 file wph 18.8.21
if($type eq "IDT"){
	my $base=judge_lastest_sitedb($types);
        my $fubase=judge_lastest_sitedb('IDTFB');
        %disease=read_disease($disease_file);
        %database=read_database($base,$type);
        %fubiaobasee=read_database_2($fubase,$type);
        &CS_get_basesite($path);
}else{
	print "$usage!"; 
}
close OUT;
close RT;
close RT2;
####################################### sub Control Deteion ####################################################
sub Deteion{
	my $shuju=shift;
	my %hash;
	open IN,"$shuju" || die "can't open the file !";
	while(<IN>){
		chomp;
		my @arr=split/\t/,$_;
		my $id=join("\t",@arr[0..3]);
		my $value="$arr[4]\t$arr[5]";
		$hash{$id}=$value;
	}
	return %hash;
}

################################## site database #############################################
sub read_database{
	my $shuju=shift;
	my $product=shift;
	open IN,"$shuju" || die "can't open the file";
	my %hash;
	my $aa=<IN>;
	chomp $aa;
	$hash{'head'}=join"\t",(split /\t/,$aa)[79..90]; #zhouwei change 20170717
    	while(<IN>){
		chomp;
		my @tmp=split/\t/,$_;
		my $key=join("_",@tmp[3,4,6,7]); #chr\tstart\t$end\tref\talt\t HGMD id as headi
		my $key2="$key\_$tmp[85]";
		my $n=@tmp;
                if($n==90){
                        $tmp[90]='-';
                }
                if($tmp[0]=~/^T1/){
                        $sub_sampleline{$key2}=1;
                }
                if($tmp[78] eq '-'){

                }
                my $info = join("\t",@tmp[79..90]);
		push @{$hash{$key}},$info ;
	}
	close IN;
	return %hash;
}

sub read_database_2{ ################### sunxq 20180408
	my $shuju=shift;
        my $product=shift;
        open IN,"$shuju" || die "can't open the file";
        my %hash;
        my $aa=<IN>;
        chomp $aa;
        $hash{'head'}=join"\t",(split /\t/,$aa)[79..90]; 
        while(<IN>){
		chomp;
                my @tmp=split/\t/,$_;
                my $key=join("_",@tmp[3,4,6,7]);
		my $info = join("\t",@tmp[78..89]);
                $info .= "\t$tmp[-3]";
                push @{$hash{$key}},$info;
                #$hash{$key}=1;
        }
        close IN;
        return %hash;
}

sub genelist{
	my $file=shift;
	my %hash;	
	my $sample;
	my @genes;#wangph change 20180606
	my $ge;#wph
	open IN,"$file" || die "can't open the file!";
	while(<IN>){
		chomp;
        	$sample=(split/\t/,$_)[0];#wph changed 20180615
        	@genes=split/, /,(split/\t/,$_)[2];#wph 180606
        	foreach $ge (@genes){
            		$hash{$sample}{$ge}=1;
		}
		#$hash{$_}=1;#wph 180606
	}
    #print Dumper(%hash);
	return %hash;
	close IN;
}

sub readsamplelist{
        my $file=shift;
        my %hash;
        open IN, "<:utf8", "$file" || die "can't open the file!";
        while(<IN>){
                chomp;
                $hash{$_}=1;
        }
        return %hash;
        close IN;
}


############################ cus filter site ###################################
sub CS_get_basesite{
	my $putin=shift;
	open IN, "$putin" || die "can't open the file!";
	my $count=1;
	while(<IN>){
		my $prefix=(split/\./,(split/\//,$_)[-1])[0]; #file list
		chomp $_;
		my $file=$_;
		open FILE,"$file" or die "can't open the file!"; # file
		chomp(my $head=<FILE>);
		if($count == 1){
			my @head_name=split/\t/,$head;
			$head_name[73]="Control_num";
			$head_name[72]="Control_ratio";
			my $names=join("\t",@head_name[0..68],"Disease system",@head_name[69..73]);
			print OUT "#Flag\tID\t$names\tnone\t$database{head}\tcase_ratio\tcase_num\tcase_ID\n";
			print RT "#Product\tID\t$names\tnone\t$database{head}\tcase_ratio\tcase_num\tcase_ID\t解读人\t末尾列\n";
			print RT2 "#Product\tID\t$names\tnone\t$database{head}\tcase_ratio\tcase_num\tcase_ID\t解读人\t末尾列\n";
		}
		while(<FILE>){
			chomp;
			my $disease_system;
			next if(/^CHROM/);
			my @tmp=split/\t/,$_;
			my @allg=split /;/,$tmp[6];
			my $exits_flag=0;
			foreach my $k(@allg){
				if(exists $gene{$prefix}{$k})#wph changed 20180615
                                {
                                        $exits_flag=1;
                                }
			}
			#next if($exits_flag==0);
			my $Ratenum=0;
			my $Ratenum05=0;
			for(19..22){
				next if($tmp[$_] eq "." || $tmp[$_]<=0.01);
				$Ratenum++;
			}
			
			for(19..22){
                                next if($tmp[$_] eq "." || $tmp[$_]<=0.05);
                                $Ratenum05++;
                        }
			
			my $ID=join("\t",@tmp[0,1,3,4]);
			if(exists $detetion{$ID}){
				my ($ratio,$num)=split/\t/,$detetion{$ID};
				$tmp[73]=$num;#control 集检出人数
				$tmp[72]=$ratio; # control 集检出率
			}else{
				$tmp[73]="-";
				$tmp[72]="-";
			}
			if(exists $disease{$tmp[6]}){
				$disease_system=$disease{$tmp[6]};
			}else{
		                $disease_system="-";
			}

			my $key2=join("_",@tmp[0,1,3,4]); #add end HGMD id 2016.08.23
        		
			my $re_info=join("\t",@tmp[0..68],$disease_system,@tmp[69..73]); 
			my $flag_2="";
			my ($flag_RT,$flag_OUT)=(0,0);
			if(exists $database{$key2}){
				$flag_2.="exists in readbase;";
			}
		
			if(exists $fubiaobasee{$key2} ){
                                $flag_2.="exists in fubiaoreadbase;";
                        }			

			
			#if((($re_info=~/trues;DM/i &&  $Ratenum05<4)|| $tmp[52] ne '-' )&& $flag_OUT==0){  
			if((($re_info=~/trues;DM/i && $Ratenum<2) || $tmp[52] ne '-' )&& $flag_OUT==0){ ### panqi 181017
				#zhowuei 表明true的DM位点被保留。
				&right_site($prefix,$re_info,$exits_flag);
				$flag_RT=1;next;
			}
			
			if($tmp[45] <4){
                                $flag_OUT=1;
                                $flag_2.="DP<4;";
			}

			if($Ratenum!=0){
				$flag_OUT=1;
				$flag_2.="AF;";
				#print OUT "AF$flag_2\t$prefix\t$re_info\n";#过滤AF
			}
			
			if($flag_OUT==0)
            		{
				&right_site($prefix,$re_info,$exits_flag);
				$flag_RT=1;next;
			}

			if($flag_2 eq ""){$flag_2="Not_in_right";}
			my $cut22=&add_sitebase_info($re_info);
        		if(!exists $case_ditacion_ratio{$key2}){$case_ditacion_ratio{$key2}="-\t-\t-";}
			print OUT "$flag_2\t$prefix\t$re_info\t$cut22\t$tmp[12]\t$case_ditacion_ratio{$key2}\n";
		}
		close FILE;
		$count++;
	}
	close IN;
	my $null="-\t" x 91;
	foreach my $key(keys %samplelist){
		if(!exists $RTsampleIDs{$key}){
			print RT "$type\t$key\t$null\n";
		}
	}
	print RT "\n\n";
	foreach my $samid(keys %dupsites){
                my $numkey=scalar (keys %{$dupsites{$samid}});
                my @key_needed = keys %{$dupsites{$samid}};
                my $ll=join(";",@key_needed);
                print RT "#$samid\t$ll\n" if($numkey>1);
        }
}

sub right_site{#add flag,RT2,add no filter gene site wph changed 2018.8.21
	my @arr=@_;
	my $prefix=$arr[0];
	my $all=$arr[1];
    	my $flag=$arr[2];
	my @tmp=split/\t/,$all;
	my $key=join("_",@tmp[0,1,3,4]);
	my $key22=join("\t",@tmp[0,1,3,4]);
	if($flag==1){
		$RTsampleIDs{$prefix}=1;
	}
	$dupsites{$key22}{$prefix}=1;
	my @allg=split /;/,$tmp[6];
	if(!exists $case_ditacion_ratio{$key22}){$case_ditacion_ratio{$key22}="-\t-\t-";}
	my %outhash;
	my $tmp3=&add_sitebase_info($all);
	if(exists $database{$key} ){  ##zhouwei  如果是false位点的话，那么解读库中的HGMD id号要写成-，否则不会进行匹配。在right文件中就会被过滤掉
		for(@{$database{$key}}){
			my @tt=split /\t/,$_;
			my $key2="$key\_$tt[6]";
			$outhash{$key2}=1;
			if(exists $sub_sampleline{$key2}){print RT "Ti";}
			my @aaa=split /\t/,$_; 
			if($aaa[11] !~ /NM/){
				$aaa[11] = $tmp[12];
				my $lll=join("\t",@aaa);
				print RT2 "Done_addNM\t$prefix\t$all\t$_\t$lll\t$case_ditacion_ratio{$key22}\n"; #add　the ID column
				if($flag == 1){
					print RT "Done_addNM\t$prefix\t$all\t$_\t$lll\t$case_ditacion_ratio{$key22}\n";
				}
			}else{
				print RT2 "Done\t$prefix\t$all\t$_\t$case_ditacion_ratio{$key22}\n"; #add　the ID column
				if($flag == 1){
					print RT "Done\t$prefix\t$all\t$_\t$case_ditacion_ratio{$key22}\n"; 	
				}
			}
		}
	}
	
	if(exists $fubiaobasee{$key}){
	    for(@{$fubiaobasee{$key}}){
		my @tt=split /\t/,$_;
            	my $key2="$key\_$tt[6]";
            	if(!exists $outhash{$key2}){
			my @aaa=split /\t/,$_; 
                        if($aaa[11] !~ /NM/){
                       		$aaa[11] = $tmp[12];
                       		my $lll=join("\t",@aaa);
				print RT2 "DoneFubiao_addNM\t$prefix\t$all\t$lll\t$case_ditacion_ratio{$key22}\n"; #add　the ID column
				if($flag == 1){
					print RT "DoneFubiao_addNM\t$prefix\t$all\t$lll\t$case_ditacion_ratio{$key22}\n";
				}
			}else{
                      		print RT2 "DoneFubiao\t$prefix\t$all\t$_\t$case_ditacion_ratio{$key22}\n"; #add　the ID column
				if($flag == 1){
					print RT "DoneFubiao\t$prefix\t$all\t$_\t$case_ditacion_ratio{$key22}\n";
				}
                    	}
                        $outhash{$key2}=1;
		}
            }
	}
	my @tmp2=split /\t/,$tmp3;
	my $key2="$key\_$tmp2[6]"; 
	if(!exists $outhash{$key2}){
		print RT2 "Todo\t$prefix\t$all\t$tmp3\t$tmp[12]\t$case_ditacion_ratio{$key22}\n";
		if($flag == 1){
			print RT "Todo\t$prefix\t$all\t$tmp3\t$tmp[12]\t$case_ditacion_ratio{$key22}\n";	
		}
		$outhash{$key2}=1;
	}
}

sub add_sitebase_info{
	my $site_info=shift;
	next if /^CHROM/;
	my @cut=split/\t/;
	my @add_info=(("-")x11);

	my $mutation_type=
            ($cut[8]=~/stopgain/i)      ? "终止"  :
            ($cut[8]=~/nonframeshift/i) ? "非移码":
            ($cut[8]=~/frameshift/i)    ? "移码"  :
            ($cut[8]=~/nonsynonymous/i) ? "错义"  :
            ($cut[8]=~/synonymous/i)    ? "同义"  : 
	    ($cut[8]=~/splicing/i)      ? "剪切"  : "?";
	
	my ($ratio_chinese,$ratio_ExAC,$ratio_gnomAD,$ratio_novo,@detection_ratio_n);
    	my  $ncon=0;
    	if($cut[20] eq "."){$ratio_chinese='未收录'; $ncon++; push @detection_ratio_n ,'千人项目东亚人群';}else{$ratio_chinese=$cut[20];}
    	if($cut[21] eq "."){$ratio_ExAC='未收录'; $ncon++;  push @detection_ratio_n, 'ExAC数据库东亚人群';}else{$ratio_ExAC=$cut[21];}
    	if($cut[19] eq "."){$ratio_gnomAD = '未收录'; $ncon++; push @detection_ratio_n, 'gnomAD数据库东亚人群';}else{$ratio_gnomAD=$cut[19];}
    	if($cut[22] eq "."){$ratio_novo='未收录';$ncon++; push @detection_ratio_n,'本地数据库';}else{$ratio_novo=$cut[22];}
    	my $ratio_detection='该变异在';
    	if($ncon==4) {
		$ratio_detection.='千人项目东亚人群、ExAC数据库东亚人群、gnomAD数据库东亚人群、本地数据库均未收录。';
	}elsif($ncon==0){
		$ratio_detection.="千人项目东亚人群检出频率为$ratio_chinese，ExAC数据库东亚人群检出频率为$ratio_ExAC，gnomAD数据库东亚人群中检出频率为$ratio_gnomAD，本地数据库检出频率为$ratio_novo。";
	}else{
		if($ratio_chinese ne '未收录'){$ratio_detection.="千人项目东亚人群检出频率为$ratio_chinese，";}
		if($ratio_ExAC ne '未收录'){$ratio_detection.="ExAC数据库东亚人群检出频率为$ratio_ExAC，";}
		if($ratio_gnomAD ne '未收录'){$ratio_detection.="gnomAD数据库东亚人群中检出频率为$ratio_gnomAD，";}
		if($ratio_novo ne '未收录'){$ratio_detection.="本地数据库检出频率为$ratio_novo，";}
		my $l=join("、",@detection_ratio_n);
		$ratio_detection.=$l."未收录。";
	}
	my $splicingds='';
	if($cut[8]=~/splicing/i){
		$splicingds="该变异经可变剪切软件进行预测，";
        if($cut[28] > 0.6 and $cut[29] > 0.6 and abs($cut[32]) >=2)
        {
            $splicingds.="DbscSNV1.1、Spidex软件预测该变异对基因或基因产物有害。";
        }
        elsif($cut[28] > 0.6 and $cut[29] > 0.6)
        {
            $splicingds.="DbscSNV1.1软件预测该变异对基因或基因产物有害。";
        }
        elsif(abs($cut[32]) >=2)
        {
            $splicingds.="Spidex软件预测该变异对基因或基因产物有害。";
        }
        else
        {
            $splicingds='';
        }
    }
=cut wph chaged 181107
		my ($ada_score,$rf_score,$dpsi_zscore);
		my $ncon1=0;
		my $ncon2=0;
		if($cut[28] eq '.' ){$ncon1++;$ada_score='无分值';}else{$ada_score=$cut[28];}
		if($cut[29] eq '.' ){$ncon1++;$rf_score='无分值';}else{$rf_score=$cut[29];}
		if($cut[32] eq '.' ){$ncon2++;$dpsi_zscore='无分值';}else{$dpsi_zscore=$cut[32];}
		if($ncon1==2 && $ncon2==1){
			$splicingds.="DbscSNV1.1和Spidex都没有预测结果。"
		}else{
			$splicingds.="DbscSNV1.1检测的ada_score分值为$ada_score,rf_score分值为$rf_score\，Spidex检测的dpsi_zscore分值为$dpsi_zscore。";
		}
	}
=cut
	my $software_pred='';
	my @softwares=split /;/,$cut[34];
	my @ssss;
	my $sss_n=0;
	if($softwares[0] =~/Deleter/ || $softwares[0] =~/dama/){push @ssss,'SIFT';$sss_n++;}  #SIFT、Polyphen2_HVAR、Polyphen2_HDIV、M-CAP
	if($softwares[1] =~/Deleter/ || $softwares[1] =~/dama/){push @ssss,'Polyphen2_HVAR';$sss_n++;}
	if($softwares[2] =~/Deleter/ || $softwares[2] =~/dama/){push @ssss,'Polyphen2_HDIV';$sss_n++;}
	if($softwares[3] =~/Deleter/ || $softwares[3] =~/dama/){push @ssss,'M-CAP';$sss_n++;}
	my $sssware=join("、",@ssss);
	if($sss_n>0){
		$software_pred=$sssware."软件预测该变异对基因或基因产物有害。"
	}
	my $wenxianzhichi='';
            #($cut[58]=~/\S/ and (length $cut[58]) >10 and $cut[59]=~/http/) ? "由于研究较少，" : "由于暂无相关研究，";
	
	my $add_info_all=();
	#my @jiyinmiaoshu=&get_jiyinmiaoshu($cut[6]);
	@add_info[1,2]=($cut[43]=~/true/i) ? ("待解读","待解读") : ("自动","临床意义未明") ; #ACMG条目,#致病性结论
	#@add_info[3,4,7,8,9]= @tmp[0,1,2,3,4];  #致病基因 疾病名称 疾病描述 遗传模式 临床指导建议 基因描述    临床表现 疾病大类
    @add_info[10,11]=@cut[13,14]; #CDS    PEP
	my $muty="$cut[14]${mutation_type}";
	$muty=~s/^-//;
    $add_info[5]="本筛查检测出$cut[6]基因的变异，变异位点为$cut[13]，查询ClinVar等公共数据库显示，$cut[13]变异会导致${muty}变异。$ratio_detection"."${software_pred}${splicingds}${wenxianzhichi}行ACMG标准，判定该变异的致病性为临床意义未明。具体情况请结合临床相关资料综合判断。";#位点描述
	#$add_info[6]="-"; #参考文献
	#push @add_info_all, join("\t", @add_info[0..11]);
	$add_info_all=join("\t", @add_info[0..11]);
	return $add_info_all;
}


sub read_disease{
	my $disease=shift;
	my %h;
	die "$disease is not exists" unless -e $disease;
	open IN,$disease or die "Can not open file $disease";
	while(<IN>){
		chomp;
		my ($gene,$system)=split/\t/,$_,2;
		$h{$gene}=$system;
	}
	close IN;
	return %h;
}

sub judge_lastest_sitedb{
	my $product_type=shift;
	my %hash_product_sitepath=();
	my @all_date;
	my $lastest_date;
	my $lastest_sitedb_path;
	my $gxy_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY"; #zhouweichange 20170712
	my $cusi_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_CuS"; #zhouweichange 20170629
        my $wes_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES";
        my $wes_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_fubiao/";
	my $genetic_sitedb_path="/ifs/TJPROJ3/HEALTH/zhongying/00research/00report/01database";
	my $cr_sitedb_path="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/CR";  #zhouweichange 20170629
	my $NNY_sitedb_path="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/Nuonanyu";  #zhouweichange 20170629
	my $genetic_black_box="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/genetic_heihe";   #zhouweichange 20170629
        my $cusi_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Fubiao_CuS/";  #sunxiaoqing 
	$hash_product_sitepath{'FWGXY'}=$gxy_sitedb_path;
	$hash_product_sitepath{'GXY'}=$gxy_sitedb_path;
	$hash_product_sitepath{'CS'}=$cusi_sitedb_path;
        $hash_product_sitepath{'IDT'}=$wes_sitedb_path;
	$hash_product_sitepath{'IDTFB'}=$wes_sitedb_path_fubiao;
	$hash_product_sitepath{'YC'}=$genetic_sitedb_path;
	$hash_product_sitepath{'YC133'}=$genetic_sitedb_path;
	$hash_product_sitepath{'CR'}=$cr_sitedb_path;
	$hash_product_sitepath{'NNY'}=$NNY_sitedb_path;
	$hash_product_sitepath{'black'}=$genetic_black_box;
        $hash_product_sitepath{'fubiao'}=$cusi_sitedb_path_fubiao;
	if( exists $hash_product_sitepath{$product_type}){
		my @all_path=glob "$hash_product_sitepath{$product_type}/*";
		foreach(@all_path){
			my$date=(split/\//)[-1];
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

sub sub_case_ditacion_ratio{
    my $detection_ratio_prefix="snp_indel_detection_out";
    my $latest_detection_file=judge_lastest_sitedb2($case_detection,$detection_ratio_prefix);#查找最新的样本检出率文件
say '*'x25;
say "正在读取最新高血压检出率文件： 
$latest_detection_file";
    open CASEDI,"$latest_detection_file" || die "can't open the $latest_detection_file!";
        %case_ditacion_ratio=map{chomp;(join"\t",(split/\t/,$_)[0..3]) => (join"\t",(split/\t/,$_)[4..6])}<CASEDI>;
    close CASEDI;
say '读取完成！！';
}

sub judge_lastest_sitedb2{
    my ($path,$prefix)=@_;
    my @all_version=();
    my $version_num;
    my $lastest_sitedb_path;
    my @all_path=glob "$path/*";
    foreach(@all_path){
        my $detection_file=basename($_);
        my $version=($detection_file=~/\_v/i) ? (split/\_/,$detection_file)[1] : "";
        if($version=~/^v(\d+)$/i){
            $version_num=$1;
            push @all_version,$version_num;
        }
    }
    (!@all_version) && (say "未发现符合条件的文件");
    my $lastest_version=max(@all_version);
    my $lastest_version_add_v="v".$lastest_version;
    my @lastest_sitedb_path=glob "$path/*\_$lastest_version_add_v\_*";
    foreach(@lastest_sitedb_path){
        if(/$prefix$/i){
            $lastest_sitedb_path=$_;
             }
    }
    return $lastest_sitedb_path;
}

