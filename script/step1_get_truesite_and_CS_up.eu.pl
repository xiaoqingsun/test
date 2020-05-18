use strict;
use warnings;
use utf8;
#use utf8::all;
use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
#use Encode;
use Spreadsheet::XLSX;
use feature qw/say/;
use File::Find;
use List::Util qw/max min sum maxstr minstr shuffle/;

#binmode(STDIN, ':encoding(utf8)');
#binmode(STDOUT, ':encoding(utf8)');
#binmode(STDERR, ':encoding(utf8)');

my ($path,$out,$types,$list);
GetOptions(
	'path=s'  => \$path,
	'out=s' => \$out,
	'type=s'	=> \$types,
	'list=s' => \$list,
);

my $usage=<<END;

perl $0 
	-path	需要检测样本的路径list文件，例如/ifs/TJPROJ3/HEALTH/Project/2C/160926_ST-E00192_0305_BH35GVALXX/cr/NHE0000000649/NHE0000000649/04YC/NHE0000000649.final.2.txt
	-type	产品类型,YC | YC133(遗传病)、GXY、CS、CR 、NNY
	-out	输出路径,输出相应产品的true的点，且YC，GXY，CS会与最新位点数据库匹配,返回与数据库匹配的结果
END

my @cut;
my %disease;
my %database;
my %oldbase;
my %fubiaobasee;
my %database_2;
my %sub_sampleline;
my %case_ditacion_ratio;
my %control_ditaction_ratio;
my %clivarPath;
my %RTsampleIDs;
my %dupsites;
my %product=(
	"CS"=>"CS",
	"GXY"=>"GXY",
);
my $type=$product{$types};
my $control_detetion="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/control_detection_ratio/170122_v0_1256wes_detection_ratio_out"; #zhouwei change 20170712 
my $cusigene="/PUBLIC/pipline/database/siteFilterDB/product/CuS/cusi_xinkang.geneV4"; #zhouwei change 20170712 
my $cusigene_info="/PUBLIC/pipline/database/siteFilterDB/product/CuS/XK_cus_disease_subtype.info"; # panqi change 20180426
my $clinvar_start_file="/PUBLIC/pipline/database/Clinvar/clivar_omim_stat.txt"; #wph update 191017
my $control_detection_path ||="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/control_detection_ratio";
my $case_detection ||="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/gxy_detection_ratio/";    #case集样本位点检出率
&sub_control_ditaction_ratio;
&sub_case_ditacion_ratio;
&read_clinvar_omim($clinvar_start_file);
my $disease_script=judge_lastest_sitedb('Disease');           ###add panqi 190709
my %gene_disease_relation=&sub_omim_id_auto($disease_script); ###add panqi 190709
my %samplelist=&readsamplelist($list);
my %gene=&genelist($cusigene);
my %detetion=&Deteion($control_detetion);
our %protein=(
	'A'=>'Ala','C'=>'Cys','D'=>'Asp','E'=>'Glu','F'=>'Phe','G'=>'Gly','H'=>'His','I'=>'Ile','K'=>'Lys','L'=>'Leu',
	'M'=>'Met','N'=>'Asn','P'=>'Pro','Q'=>'Gln','R'=>'Arg','S'=>'Ser','T'=>'Thr','V'=>'Val','W'=>'Trp','Y'=>'Tyr',
	'X'=>'*',
);  #wph add 190122

################################################# main #######################################################################
open OUT, ">$out/$types.err.txt" || die "can't open the $out/$types.err.txt!";
open RT, "> $out/$types.right.txt" || die "can't open the $out/$types.right.txt!";
binmode(RT, ':encoding(utf8)');
binmode(OUT, ':encoding(utf8)');
if($type=~/GXY/){
	my $base=judge_lastest_sitedb($type);
	%database=read_database($base,$type);
	%disease=read_disease($cusigene_info);
	&get_basesite($type,$path);
}elsif($type eq "CS"){
	my $base=judge_lastest_sitedb($type);
	my $fubase=judge_lastest_sitedb('CSFB');
	%disease=read_disease($cusigene_info);
	%database=read_database($base,$type);
	%oldbase=&read_olddatabase;
	%fubiaobasee=read_database_2($fubase,$type);
	&CS_get_basesite($path);
}else{
	print "$usage!"; 
}
close OUT;
close RT;

####################################### sub Control Deteion ####################################################
sub Deteion{
	my $shuju=shift;
	my %hash;
	open IN, "<:utf8", "$shuju" || die "can't open the file !";
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
	open IN, "<:utf8", "$shuju" || die "can't open the file";
	my %hash;
	my $aa=<IN>;
	chomp $aa;
	my @tmp=split/\t/,$aa;
	$hash{'head'}=(join"\t",@tmp[79..90]); #zhouwei change 20170717
    	while(<IN>){
		chomp;
		my @tmp=split/\t/,$_;
		my $key=join("_",@tmp[3,4,6,7]); 
		my $key2="$key\_$tmp[85]"; 
		my $n=@tmp;
		if($n==90){
			$tmp[90]='-';
		}
		if($tmp[0]=~/^T1/){
			$sub_sampleline{$key2}=1;
		}
		
		if($product=~/GXY/){
			my $info = join("\t",@tmp[45,67..78]);
			push @{$hash{$key}},$info;
		}else{
			my $info = join("\t",@tmp[79..90]);
			push @{$hash{$key}},$info ;
		}
	}
	close IN;
	return %hash;
}

sub read_database_2{ ################### sunxq 20180408
	my $shuju=shift;
	my $product=shift;
	open IN, "<:utf8", "$shuju" || die "can't open the file";
	my %hash;
	#my $key;
	#my @tmp=();
	my $aa=<IN>;
	chomp $aa;
	$hash{'head'}=join"\t",(split /\t/,$aa)[79..90];
	while(<IN>){
		chomp;
		my @tmp=split/\t/,$_;
		my $key=join("_",@tmp[3,4,6,7]);
		my $info=join("\t",@tmp[79..90]);
		#$info .= "\t$tmp[-3]"; #wph del 190327
		push @{$hash{$key}},$info;
	}
	close IN;
	return %hash;
}


sub read_olddatabase{
	my $shuju="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_CuS/201711221046/XK_cusbase201711221046.right.txt";
	my %hash_o;
	open IN,"<:utf8", "$shuju" || die "can't open the file";
	while(<IN>){
		chomp;
		my @tmp=split/\t/,$_;
		my $key=join("_",@tmp[4,5,7,8]); 
		my $info=join("\t",@tmp[67,68,71,72]); 
		$hash_o{$key}=$info;
	}
	close IN;
	return %hash_o;
}

sub genelist{
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
		open FILE,"<:utf8","$file" or die "can't open the file!"; # file wangpenghui 181017
		chomp(my $head=<FILE>);
		if($count == 1){
			$head =~s/True_Flase_HGMD/HGMD_True_Flase/g;
			$head =~s/True_Flase_Clinvar/Clinvar_True_Flase/g;
			$head =~s/true_genome_hgmd_contex/hgmd_true_genome_contex/g;
			$head =~s/true_pep_hgmd_contex/hgmd_true_pep_contex/g;
			$head =~s/true_genome_clivar_contex/clivar_true_genome_contex/g;
			$head =~s/true_pep_clivar_contex/clivar_true_pep_contex/g;#wph chang head 190329
			my @head_name=split/\t/,$head;
			$head_name[73]="Control_num";
			$head_name[72]="Control_ratio";
			my $names=join("\t",@head_name[0..68],"Disease class",@head_name[69..73]);
			print OUT "#Flag\tID\t$names\tolddatabase\t$database{head}\tcase_ratio\tcase_num\tcase_ID\n";
			print RT "#Product\tID\t$names\tolddatabase\t$database{head}\tcase_ratio\tcase_num\tcase_ID\t解读人\t末尾列\n";
		}
		while(<FILE>){
			chomp;
			my $disease_system='-';
			next if(/^CHROM/);
			my @tmp=split/\t/,$_;
			my @allg=split /;/,$tmp[6];
			my $exits_flag=0;
			foreach my $k(@allg){
				if(exists $gene{$k})
				{
					$exits_flag=1;
				}
			}
			next if($exits_flag==0);
			my $Ratenum=0;
			my $Ratenum05=0;
			for(19..22){
				next if($tmp[$_] eq "." || $tmp[$_]<=0.001);
				$Ratenum++ if($_ != 20);
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
				$tmp[73]="0"; #wph chaged 190103
				$tmp[72]="0";
			}
			
			my $key2=join("_",@tmp[0,1,3,4]); #add end HGMD id 2016.08.23
        		
			my $re_info=join("\t",@tmp[0..68],$disease_system,@tmp[69..73]); 
			my $flag_2="";
			my $flag_22="";
			my ($flag_RT,$flag_OUT,$flag_DP)=(0,0,0);
			if(exists $database{$key2}){
			#	$flag_2.="exists_in_readbase;";
			}
			if(exists $fubiaobasee{$key2} ){
			#	$flag_2.="exists_in_fubiaoreadbase;";
			}	
		
			if($tmp[45] <10){
				$flag_DP=1;
				$flag_2.="DP<10;";
				
				#print OUT "DP$flag_2\t$prefix\t$re_info\n";#所有点过DP、AF
			}
			if($tmp[52] =~/negative/){
				$flag_DP=1;
				$flag_2.="Del_site_hotspot;";
				
			}
			
			#if((($re_info=~/trues;DM/i && $Ratenum05<4) || $tmp[52] ne '-' )&& $flag_OUT==0){  
			if( $tmp[52] ne '-' && $flag_DP==0 && $tmp[52] !~/negative/){ ### panqi 181017,181109
				#zhowuei 表明true的DM位点被保留。
				&right_site($prefix,$re_info,0,$flag_2,0);
				$flag_RT=1;next;
			} ### remove ($re_info=~/trues;DM/i && $Ratenum<2) panqi 190124

			if($Ratenum!=0){
				$flag_OUT=1;
				$flag_2.="AF;";
			
				#print OUT "AF$flag_2\t$prefix\t$re_info\n";#过滤AF
			}
			
			my @protoin_predict=split/;/,$tmp[34];
			my ($SIFT_Polyphen2,$CADD)=(@protoin_predict>=4) ? ( (join"",@protoin_predict[0..3]),$protoin_predict[-1] ) : ($tmp[34],$tmp[34]); ##changed by sunxq 20171122
			#if(@protoin_predict>=5 and $SIFT_Polyphen2!~/deleterious|damaging/i and $SIFT_Polyphen2 ne "----" )
			#{ $flag_OUT=1;$flag_2.="protoin_predict_filter;";}
			
			my $protoin_predict_num=0;
			for(@protoin_predict){
				if($_=~/tolerated|benign/i){$protoin_predict_num++;}
				if($_=~/deleterious|damaging/i){$protoin_predict_num=-1;} #wph add 190115
			}
			if($protoin_predict_num>=3){$flag_OUT=1; $flag_2.="protion_predict_filter;";}

			my $flag_dm;
			my $flag_dm_type;
			my ($FLAG_DM_62,$FLAG_DM_64,$FLAG_DM_T62,$FLAG_DM_T64)=(0,0,0,0);
			my @flag_dm_array;
			for(my $i=62;$i<66;$i++)
			{
				next if($i==63 || $i==65);
				next if($tmp[$i] eq "-");
				@flag_dm_array=split /\|/,$tmp[$i];
				for(@flag_dm_array)
				{
					$flag_dm_type=(split /:/,$_)[0];
					next if($flag_dm_type!~/true_type/);
					$FLAG_DM_T62=1 if($i==62);
					$FLAG_DM_T64=1 if($i==64);

					$flag_dm=(split /:/,$_)[3];
					if($flag_dm!~/DM|pathogenic|significance/i && $flag_dm ne "not provided" && $flag_dm ne "other" )
					{
						$FLAG_DM_62=2 if($i==62);
						$FLAG_DM_64=2 if($i==64);
					}
					
					if($flag_dm=~/DM|pathogenic/i  && $flag_dm!~/conflict/i )
					{
						$FLAG_DM_62=1 if($i==62);
						$FLAG_DM_64=1 if($i==64);
					}
				}
			}
			
			my $mu_type=0;
			my @acmg=&acmg_auto_judge($re_info,$prefix);
			if($tmp[8] eq "synonymous SNV" || ($tmp[5]!~/splicing/ && $tmp[5]!~/exonic/) )
			{
				$mu_type=1;
				if( ($FLAG_DM_T62==1 && $FLAG_DM_62==1) || ( $FLAG_DM_T64==1 && $FLAG_DM_64==1) )
				{
					$flag_RT=1; $flag_22.="DM;"
				}
			}
			else{
				if(@acmg>0){
					$flag_RT=1;
					for(@acmg){$flag_22.="$_;"}
				}
			}
			
			#if(($tmp[5] eq 'intronic' || $tmp[8] eq 'synonymous SNV' ) && $FLAG_DM_T62==0 && $FLAG_DM_T64==0 ){ #sxq changed 190508
			#	$flag_OUT=1;
			#	$flag_2.="intronic_syn;";
			#}
			
			if($flag_DP==0){
				if( ($flag_RT==1 && $mu_type==1 && $flag_OUT==0)  || ($flag_RT==1 && $mu_type==0) || ($mu_type==0 && $flag_OUT==0) )
				{
					if($flag_RT==1 ){
						&right_site($prefix,$re_info,0,$flag_22,$FLAG_DM_64); 
					}
					elsif($flag_OUT==0){
						&right_site($prefix,$re_info,0,$flag_2,$FLAG_DM_64);
					}
				}
				else{
					&right_site($prefix,$re_info,1,$flag_2,$FLAG_DM_64);
				}
			}
			elsif($flag_DP==1){
				&right_site($prefix,$re_info,1,$flag_2,$FLAG_DM_64);
			}
			
			
=head			
			if( ($flag_RT==1 && $flag_22 ne "") || ($flag_RT==0 && $flag_2 eq "") )
            		{
				&right_site($prefix,$re_info,0,$flag_2);
			}
			else{
				if($flag_RT==1){
					&right_site($prefix,$re_info,1,$flag_22);
				}
				elsif($flag_RT==0){
					&right_site($prefix,$re_info,1,$flag_2);
				}
			}

			if($flag_2 eq ""){$flag_2="Not_in_right";}
			($re_info,my @cut)=&add_sitebase_info($re_info);
			for(@cut){
				my @arry=split /\t/,$_;
				my $diseasename=$arry[6];
				my $diseaseclass='-';
				foreach my $k(@allg){
					if(exists $disease{$diseasename}{$k}){
						if($diseaseclass eq '-'){$diseaseclass='';}
						if($diseaseclass ne '-' && $diseaseclass ne ''){$diseaseclass.=';';}
						$diseaseclass.="$disease{$diseasename}{$k}";	
					}
				}
				my $key22=join("\t",@tmp[0,1,3,4]); #wph add 20181112
				my $re_info2=join("\t",@tmp[0..68],$diseaseclass,@tmp[69..73]);
        			if(!exists $case_ditacion_ratio{$key22}){$case_ditacion_ratio{$key22}="0\t0\t-";}#wph chaged 190103
				print OUT "$flag_2\t$prefix\t$re_info\t-\t$_\t$tmp[12]\t$case_ditacion_ratio{$key22}\n"; #del ref_info2
			}
=cut
			
		}
		close FILE;
		$count++;
	}
	#close OUT;
	#close RT;
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

sub get_jiyinmiaoshu{
	my $gene=shift;
	my $cha_gene=$gene=~s/\(.*\)//gr;
	my @all_gene=split/;/,$cha_gene;
	my %hash=();
	#my @cut;
	my @all_term;

	open INN,"<:utf8", $cusigene_info or die $!;
	##<INN>;
	while(<INN>)
	{
		chomp;
		@cut=split /\t/,$_;
		push @{$hash{$cut[1]}}, join("\t",@cut[7,6,3,4,1]);
	}
	close INN;

	foreach my $eachgene (@all_gene){
		foreach my $item ( @{$hash{$eachgene}} ){
			push @all_term, $item;
		}
	}
	
	return @all_term;
}

sub add_sitebase_info{
	my $site_info=shift;
	next if /^CHROM/;
	my @cut=split/\t/,$site_info;
	my @add_info=(("-")x11);

	my $extra_type="";
	my $mutation_type=
            ($cut[8]=~/stopgain/i)      ? "无义"  :
	    ($cut[8]=~/stoploss/i)	? "终止缺失":
            ($cut[8]=~/nonframeshift/i) ? "非移码":
            ($cut[8]=~/frameshift/i)    ? "移码"  :
            ($cut[8]=~/nonsynonymous/i) ? "错义"  :
            ($cut[8]=~/synonymous/i)    ? "同义"  : 
	    ($cut[8]=~/splicing/i)    ? "剪切"  : "?";
	
	my ($ratio_chinese,$ratio_ExAC,$ratio_gnomAD,$ratio_novo,@detection_ratio_n);
    my  $ncon=0;
    if($cut[20] eq "."){$ratio_chinese='未收录'; $ncon++; push @detection_ratio_n ,'千人项目全球人群';}else{$ratio_chinese=$cut[20];}
    if($cut[21] eq "."){$ratio_ExAC='未收录'; $ncon++;  push @detection_ratio_n, 'ExAC数据库欧洲人群';}else{$ratio_ExAC=$cut[21];}
    if($cut[19] eq "."){$ratio_gnomAD = '未收录'; $ncon++; push @detection_ratio_n, 'gnomAD数据库欧洲人群';}else{$ratio_gnomAD=$cut[19];}
    if($cut[22] eq "."){$ratio_novo='未收录';$ncon++; push @detection_ratio_n,'本地数据库';}else{$ratio_novo=$cut[22];}
    my $ratio_detection='该变异在';
    if($ncon==4) {
		$ratio_detection.='千人项目欧洲人群、ExAC数据库欧洲人群、gnomAD数据库欧洲人群、本地数据库均未收录。';
	}elsif($ncon==0){
		$ratio_detection.="千人项目全球人群检出频率为$ratio_chinese，ExAC数据库欧洲人群检出频率为$ratio_ExAC，gnomAD数据库欧洲人群中检出频率为$ratio_gnomAD，本地数据库检出频率为$ratio_novo。";
	}else{
		if($ratio_chinese ne '未收录'){$ratio_detection.="千人项目全球人群检出频率为$ratio_chinese，";}
		if($ratio_ExAC ne '未收录'){$ratio_detection.="ExAC数据库欧洲人群检出频率为$ratio_ExAC，";}
		if($ratio_gnomAD ne '未收录'){$ratio_detection.="gnomAD数据库欧洲人群中检出频率为$ratio_gnomAD，";}
		if($ratio_novo ne '未收录'){$ratio_detection.="本地数据库检出频率为$ratio_novo，";}
		my $l=join("、",@detection_ratio_n);
		$ratio_detection.=$l."未收录。";
	}

	
	
	my $splicingds='';
	if( ($cut[5]=~/splicing/i && $cut[5]!~/exon/) || $cut[8]=~/splicing/  || $cut[8]=~/^synonymous SNV/i ){  ### panqi 190714
		$splicingds="该变异经可变剪切软件进行预测，";
		my ($ada_score,$rf_score,$dpsi_zscore);
        if($cut[28] ne "." && $cut[29] ne "." && $cut[32] ne "." && $cut[28] > 0.6 && $cut[29] > 0.6 && abs($cut[32]) >=2)
        {
            $splicingds.="dbscSNV1.1、SPIDEX软件预测该变异对基因或基因产物有害。";
        }
        elsif( $cut[28] ne "." and $cut[29] ne "." and $cut[28] > 0.6 and $cut[29] > 0.6)
        {
            $splicingds.="dbscSNV1.1软件预测该变异对基因或基因产物有害。";
        }
        elsif( $cut[32] ne "." and abs($cut[32])>=2 )
        {
            $splicingds.="SPIDEX软件预测该变异对基因或基因产物有害。";
        }
        else
        {
            $splicingds='';
        }
    }
=cut wph changed 181107
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
=cut
	
	my $software_pred='';
	my @softwares=split /;/,$cut[34];
	my (@ssss,@sssk,@ssst);
	my ($sss_n,$ssk_n,$sst_n)=(0,0,0);
	if($softwares[0] =~/Deleter/i || $softwares[0] =~/dama/i){push @ssss,'SIFT';$sss_n++;}  #SIFT、Polyphen2_HVAR、Polyphen2_HDIV、M-CAP
	elsif($softwares[0] eq "-"){push @sssk,'SIFT';$ssk_n++;}
	elsif($softwares[0] !~/Deleter/i && $softwares[0] !~/dama/i){push @ssst,'SIFT';$sst_n++;}
	
	if($softwares[1] =~/Deleter/i || $softwares[1] =~/dama/i){push @ssss,'Polyphen2_HVAR';$sss_n++;}
	elsif($softwares[1] eq "-"){push @sssk,'Polyphen2_HVAR';$ssk_n++;}
	elsif($softwares[1] !~/Deleter/i && $softwares[1] !~/dama/i){push @ssst,'Polyphen2_HVAR';$sst_n++;}
	
	if($softwares[2] =~/Deleter/i || $softwares[2] =~/dama/i){push @ssss,'Polyphen2_HDIV';$sss_n++;}
	elsif($softwares[2] eq "-"){push @sssk,'Polyphen2_HDIV';$ssk_n++;}
	elsif($softwares[2] !~/Deleter/i && $softwares[2] !~/dama/i){push @ssst,'Polyphen2_HDIV';$sst_n++;}
	
	if($softwares[3] =~/Deleter/i || $softwares[3] =~/dama/i){push @ssss,'M-CAP';$sss_n++;}
	elsif($softwares[3] eq "-"){push @sssk,'M-CAP';$ssk_n++;}
	elsif($softwares[3] !~/Deleter/i && $softwares[3] !~/dama/i){push @ssst,'M-CAP';$sst_n++;}
	
	my $sssware=join("、",@ssss);
	if($sss_n>0){
		$software_pred.=$sssware."软件预测该变异对基因或基因产物有害。";
	}
	$sssware=join("、",@ssst);
	if($sst_n>0){
#		$software_pred.=$sssware."软件预测该变异对基因或基因产物无害。";
	}
	$sssware=join("、",@sssk);
	if($ssk_n>0 && $ssk_n<4){
#		$software_pred.='该变异在蛋白功能预测软件'.$sssware."中未收录。";
	}
	if($ssk_n==4 && $cut[8]=~/^nonsynonymous SNV/i)
	{
#		$software_pred.='该变异在蛋白功能预测软件'.$sssware."中未收录。";
	}
	
	my $wenxianzhichi='';
            #($cut[58]=~/\S/ and (length $cut[58]) >10 and $cut[59]=~/http/) ? "由于研究较少，" : "由于暂无相关研究，";
	
	my @add_info_all=();
	my @jiyinmiaoshu=&get_jiyinmiaoshu($cut[6]);
	
	for(my $i=0;$i<@jiyinmiaoshu;$i++)
	{ 
		#my @tmp=split /\t/, encode('utf8', $jiyinmiaoshu[$i]);
		#my @tmp=split /\t/, Encode::_utf8_on($jiyinmiaoshu[$i]);
		#my $flag = utf8::is_utf8($tmp[0]);
		
		my $arhet="";
		my @tmp=split /\t/, $jiyinmiaoshu[$i];
		@add_info[0,1]=($cut[43]=~/true/i) ? ("待解读","待解读") : ("自动","临床意义未明") ; #ACMG条目,#致病性结论
		@add_info[2,3,6,7,8]= @tmp[0,1,2,3,4];  #致病基因 疾病名称 疾病描述 遗传模式 临床指导建议 基因描述    临床表现 疾病大类
		@add_info[9,10]=@cut[13,14];
		($add_info[10],$extra_type)=pep_trans($add_info[10]); #CDS    PEP
		if($cut[30]=~/het/ and $tmp[3]=~/AR/){$arhet="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，";}
		my $muty="";#wph changed 1910104 $cut[14]${mutation_type}";
		if($mutation_type eq "?" or $cut[8] =~/unknown/i){  #wph add 190418
			$mutation_type=$extra_type;
			$cut[8] =
			($extra_type eq "错义")					? "nonsynonymous SNV"	    :
			($extra_type eq "同义")					? "synonymous SNV"	    :
			($extra_type eq "无义")					? "stopgain"		    :
			($extra_type eq "终止缺失")				? "stoploss"		    :
($extra_type eq "剪切" && $cut[13] !~/del|->|>-/i && $cut[5] ne "intronic")	? "splicing SNV"	    :
			($extra_type eq "剪切" && $cut[13] =~/del|->|>-/i)	? "splicing INDEL"	    :
			($extra_type eq "非移码" && $cut[13] =~/del/) 		? "nonframeshift deletion"  :
			($extra_type eq "非移码" && $cut[13] =~/dup/) 		? "nonframeshift insertion" :
			($extra_type eq "移码" && $cut[13] =~/del/)		? "frameshift deletion"	    :
			($extra_type eq "移码" && $cut[13] =~/dup/)		? "frameshift insertion"    : "SNV";

		}
		if ($cut[5]=~/splicing/ or $extra_type eq "剪切"){
			if($cut[5] =~ /UTR/) { #wph add 190402
				$muty="位于UTR区";
			}elsif($cut[5] =~ /intronic/) {
				$muty="位于内含子区";
			}else{
				$muty="位于剪切区";
			}
		}else{
			$muty="会导致$add_info[10]${mutation_type}变异";
		}
		$muty=~s/^-//;
		$add_info[4]="本筛查检测出$cut[6]基因的变异，变异位点为$cut[13]，查询ClinVar等公共数据库显示，$cut[13]变异${muty}。$ratio_detection"."${software_pred}${splicingds}${wenxianzhichi}行ACMG标准，判定该变异的致病性为临床意义未明。${arhet}具体情况请结合临床相关资料综合判断。";#位点描述
		$add_info[5]="-"; #参考文献
		
		push @add_info_all, join("\t", @add_info[0..10]);
	}
	$site_info=join("\t",@cut);
	return $site_info,@add_info_all;
}

sub pep_trans{ #wph add 190122 190402 add type
	my $p=shift;
	my $type2="";
	if($p=~/p\.([A-Y]+)(\d+)([A-Y]+)(.*)/g){
		my $locaa=$2;
		my $firstaa =join"",(map{$protein{$_}//'-'}(split//,$1));
		my $secondaa=join"",(map{$protein{$_}//'-'}(split//,$3));
		my $lastaa=$4;
		if($p =~/fs/){
			$type2="移码";
		}elsif($firstaa eq $secondaa){
			$type2="同义";
		}elsif($secondaa =~ /\*/){
			$type2="无义";
		}elsif($firstaa =~ /\*/){
			$type2="终止缺失";
		}else{
			$type2="错义";
		}
		$p="p\.${firstaa}${locaa}${secondaa}${lastaa}";
	}elsif($p=~/p\.([A-Y]+)(\d+)fs/g){
		$type2="移码";
		my $locaa=$2;
		my $firstaa =join"",(map{$protein{$_}//'-'}(split//,$1));
		$p="p\.${firstaa}${locaa}fs";
	}elsif($p eq "-"){
		$type2="剪切";
	}elsif($p =~/fs/i){
		$type2="移码";
	}elsif($p =~/X/i){
		$type2="无义";
	}else{
		$type2="非移码";
		#return $p,$type;       
	}
	return $p,$type2;
}

sub right_site{
	my @arr=@_;
	my $prefix=$arr[0];
	my $all=$arr[1];
	my $flag_RT=$arr[2];
	my $flag_wrong=$arr[3];
	my $flag_DM=$arr[4];
	
	($all,my @tmp3)=&add_sitebase_info($all);
	my @tmp=split/\t/,$all;
	my $key=join("_",@tmp[0,1,3,4]);
	my $key22=join("\t",@tmp[0,1,3,4]);

	if( $tmp[5]=~/intronic/ || $tmp[5]=~/UTR/ || $tmp[8] eq 'synonymous SNV' )
	{
		$flag_wrong.="intronic;" if($tmp[5]=~/intronic/);
		$flag_wrong.="UTR;" if($tmp[5]=~/UTR/ );
		$flag_wrong.="synonymous_SNV;" if( $tmp[8] eq 'synonymous SNV' );
	}
	$flag_wrong.="Clinvar_Benign" if($flag_DM==2);

	$RTsampleIDs{$prefix}=1 if($flag_RT==0);
	$dupsites{$key22}{$prefix}=1 if($flag_RT==0);
	my @allg=split /;/,$tmp[6];

	if(!exists $case_ditacion_ratio{$key22}){$case_ditacion_ratio{$key22}="0\t0\t-";} #wph chaged 190103
	my %outhash;
	#($all,my @tmp3)=&add_sitebase_info($all); wph changed 190424  修改突变类型，先赋值
	if(exists $database{$key} ){  ##zhouwei  如果是false位点的话，那么解读库中的HGMD id号要写成-，否则不会进行匹配。在right文件中就会被过滤掉
		for(@{$database{$key}}){
			my @tt=split /\t/,$_;
			my $line=encode('utf8',$_);
			my @ttt=split /\t/,$line;
			my $key2="$key\_$tt[6]";
			next if(exists $outhash{$key2} );
			$outhash{$key2} = 1;
			my @anno_info=split /\t/,$all;
			print "site_info:$anno_info[30]\t$tt[5]\t$tt[8]\n";#wph changed 180104 dababase_site
			if(exists $sub_sampleline{$key2}){print RT "Ti";}
			my $diseasename=$tt[6];
			my $diseaseclass='-';
			foreach my $k(@allg){
				if(exists $disease{$diseasename}{$k}){
					if($diseaseclass eq '-'){$diseaseclass='';}
					if($diseaseclass ne '-' && $diseaseclass ne ''){$diseaseclass.=';';}
					$diseaseclass.="$disease{$diseasename}{$k}";    
				}
				if($tmp[70]=~/^chr/ && exists $gene_disease_relation{$k}{$ttt[6]} ){  ###add panqi 190709
					$tmp[70]=$gene_disease_relation{$k}{$ttt[6]};
				}
				
			}
			$tmp[69]=$diseaseclass;

			my $all2=join("\t",@tmp);
			my @aaa=split /\t/,$_; 
			if($aaa[11] !~ /NM/){
				$aaa[11] = $tmp[12];
				my $lll=join("\t",@aaa);
				print RT "Done_addNM;"."$flag_wrong\t$prefix\t$all2\t-\t$lll\t$case_ditacion_ratio{$key22}\n" if($flag_RT==0); #add　the ID column
				print OUT "$flag_wrong;exists_in_readbase;\t$prefix\t$all2\t-\t$lll\t$case_ditacion_ratio{$key22}\n" if($flag_RT==1);
			}else{
				print RT "Done;"."$flag_wrong\t$prefix\t$all2\t-\t$_\t$case_ditacion_ratio{$key22}\n" if($flag_RT==0); #add　the ID column
				print OUT "$flag_wrong;exists_in_readbase;\t$prefix\t$all2\t-\t$_\t$case_ditacion_ratio{$key22}\n" if($flag_RT==1);
			}
		}
	}
    if(exists $fubiaobasee{$key}){
		for(@{$fubiaobasee{$key}}){
			my @tt=split /\t/,$_;
			my $line=encode('utf8',$_);
			my @ttt=split /\t/,$line;
			my $key2="$key\_$tt[6]";
			if(!exists $outhash{$key2}){
				my $diseasename=$tt[6];
				my $diseaseclass='-';
				foreach my $k(@allg){
					if(exists $disease{$diseasename}{$k}){
						if($diseaseclass eq '-'){$diseaseclass='';}
						if($diseaseclass ne '-' && $diseaseclass ne ''){$diseaseclass.=';';}
						$diseaseclass.="$disease{$diseasename}{$k}";    
					}
					if($tmp[70]=~/^chr/ && exists $gene_disease_relation{$k}{$ttt[6]} ){  ###add panqi 190709
						$tmp[70]=$gene_disease_relation{$k}{$ttt[6]};
					}					
				}
				$tmp[69]=$diseaseclass;

				my $all2=join("\t",@tmp);
				my @aaa=split /\t/,$_; 
				if($aaa[11] !~ /NM/){
					$aaa[11] = $tmp[12];
					my $lll=join("\t",@aaa);
					print RT "DoneFubiao_addNM;"."$flag_wrong\t$prefix\t$all2\t-\t$lll\t$case_ditacion_ratio{$key22}\n" if($flag_RT==0); #add　the ID column
					print OUT "$flag_wrong;exists_in_fubiaoreadbase;\t$prefix\t$all2\t-\t$lll\t$case_ditacion_ratio{$key22}\n" if($flag_RT==1);
				}else{
					print RT "DoneFubiao;"."$flag_wrong\t$prefix\t$all2\t-\t$_\t$case_ditacion_ratio{$key22}\n" if($flag_RT==0); #add　the ID column
					print OUT "$flag_wrong;exists_in_fubiaoreadbase;\t$prefix\t$all2\t-\t$_\t$case_ditacion_ratio{$key22}\n" if($flag_RT==1);
				}
				$outhash{$key2}=1;
			}
		}
    }
	if(exists $oldbase{$key} ){
		my ($a1,$a2,$a3,$a4)=(split /\t/, $oldbase{$key})[0,1,2,3];
		for(@tmp3)
		{
			my @tmp2=split /\t/,$_;
			my $line=encode('utf8',$_);
			my @tmp22=split /\t/,$line;
			my $key2="$key\_$tmp2[6]"; 
			if(!exists $outhash{$key2}){
				my $diseasename=$tmp2[6];
				my $diseaseclass='-';
				foreach my $k(@allg){
					if(exists $disease{$diseasename}{$k}){
						if($diseaseclass eq '-'){$diseaseclass='';}
						if($diseaseclass ne '-' && $diseaseclass ne ''){$diseaseclass.=';';}
						$diseaseclass.="$disease{$diseasename}{$k}";
					}
					if($tmp[70]=~/^chr/ && exists $gene_disease_relation{$k}{$tmp22[6]} ){  ###add panqi 190709
						$tmp[70]=$gene_disease_relation{$k}{$tmp22[6]};
					}					
					
				}
				$tmp[69]=$diseaseclass;
				
				my $all2=join("\t",@tmp);
				print RT  "Todo_oldBD;"."$flag_wrong\t$prefix\t$all2\t$a3\t$a1\t$a2\t$tmp2[2]\t$tmp2[3]\t$tmp2[4]\t$a4\t$tmp2[6]\t$tmp2[7]\t$tmp2[8]\t$tmp2[9]\t$tmp2[10]\t$tmp[12]\t$case_ditacion_ratio{$key22}\n" if($flag_RT==0);
				print OUT "$flag_wrong\t$prefix\t$all2\t$a3\t$a1\t$a2\t$tmp2[2]\t$tmp2[3]\t$tmp2[4]\t$a4\t$tmp2[6]\t$tmp2[7]\t$tmp2[8]\t$tmp2[9]\t$tmp2[10]\t$tmp[12]\t$case_ditacion_ratio{$key22}\n" if($flag_RT==1);
				$outhash{$key2}=1;
			}
		}
	}
	for(@tmp3)
	{
		#Encode::_utf8_on($_);
		my @tmp2=split /\t/,$_;
		my $line=encode('utf8',$_);
		my @tmp22=split /\t/,$line;
		my $key2="$key\_$tmp2[6]"; 
		if(!exists $outhash{$key2}){
			my $diseasename=$tmp2[6];
			my $diseaseclass='-';
			foreach my $k(@allg){
				if(exists $disease{$diseasename}{$k}){
					if($diseaseclass eq '-'){$diseaseclass='';}
					if($diseaseclass ne '-' && $diseaseclass ne ''){$diseaseclass.=';';}
					$diseaseclass.="$disease{$diseasename}{$k}";
				}
				if($tmp[70]=~/^chr/ && exists $gene_disease_relation{$k}{$tmp22[6]} ){  ###add panqi 190709
					$tmp[70]=$gene_disease_relation{$k}{$tmp22[6]};
				}				
			}
			$tmp[69]=$diseaseclass;
			
			my $all2=join("\t",@tmp);
			print RT "Todo;"."$flag_wrong\t$prefix\t$all2\t-\t$_\t$tmp[12]\t$case_ditacion_ratio{$key22}\n" if($flag_RT==0);
			print OUT "$flag_wrong\t$prefix\t$all2\t-\t$_\t$tmp[12]\t$case_ditacion_ratio{$key22}\n" if($flag_RT==1);
			$outhash{$key2}=1;
		}
	}
}

sub sub_omim_id_auto{                                        ###add panqi 190709
	my $disease_script=shift;
	my %gene_disease_relation=();
	my @cut;
	my @arr;
	my $disease_name;
	my $tmp2;

	open INN, $disease_script or die $!;
	while(<INN>){
		chomp;
		next if($_=~/^#NewID/);
		#my $line=encode('utf8',$_);
		@cut=split /\t/,$_;
		if($cut[0] eq "OMIM:185500")
		{
			$cut[0]="OMIM:185500";
		}
		$disease_name= (split /【/,$cut[4])[0]; 		
		$gene_disease_relation{$cut[3]}{$disease_name}=$cut[0];
		@arr=split /\|/,$cut[12] if($cut[12] ne "-" && $cut[12] ne "NULL" && $cut[12] ne "");
		next if($cut[12] eq "-" || $cut[12] eq "NULL" || $cut[12] eq "");
		for(@arr)
		{
			$disease_name=(split /\_/, $_)[1];
			$gene_disease_relation{$cut[3]}{$disease_name}=$cut[0];
		}
	}
	close INN;
	return %gene_disease_relation;
}

sub read_disease{
	my $disease=shift;
	my %h;
	die "$disease is not exists" unless -e $disease;
	open IN,"<:utf8",$disease or die "Can not open file $disease";
	while(<IN>){
		chomp;
		my @aa=split/\t/,$_;
		$h{$aa[3]}{$aa[1]}=$aa[2];
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
	my $genetic_sitedb_path="/ifs/TJPROJ3/HEALTH/zhongying/00research/00report/01database";
	my $cr_sitedb_path="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/CR";  #zhouweichange 20170629
	my $NNY_sitedb_path="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/Nuonanyu";  #zhouweichange 20170629
	my $genetic_black_box="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/genetic_heihe";   #zhouweichange 20170629
	my $cusi_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_CuS_fubiao";  #sunxiaoqing
	my $gene_descript_path="/PUBLIC/pipline/database/Auto_pipline/gene_description/";#wangpenghui        ###add panqi 190709
	my $disease_descript_path="/PUBLIC/pipline/database/Auto_pipline/disease_description/";#wangpenghui  ###add panqi 190709
	$hash_product_sitepath{'FWGXY'}=$gxy_sitedb_path;
	$hash_product_sitepath{'GXY'}=$gxy_sitedb_path;
	$hash_product_sitepath{'CS'}=$cusi_sitedb_path;
	$hash_product_sitepath{'YC'}=$genetic_sitedb_path;
	$hash_product_sitepath{'YC133'}=$genetic_sitedb_path;
	$hash_product_sitepath{'CR'}=$cr_sitedb_path;
	$hash_product_sitepath{'NNY'}=$NNY_sitedb_path;
	$hash_product_sitepath{'black'}=$genetic_black_box;
	$hash_product_sitepath{'CSFB'}=$cusi_sitedb_path_fubiao;
	$hash_product_sitepath{'Gene'}=$gene_descript_path;               ###add panqi 190709
	$hash_product_sitepath{'Disease'}=$disease_descript_path;         ###add panqi 190709
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
###################################### 阴性对照检出率 #################################################
sub sub_control_ditaction_ratio{
    my $control_detection_prefix="detection_ratio_out";
    my $latest_control_detection_ratio=judge_lastest_sitedb2($control_detection_path,$control_detection_prefix);
say '*'x25;
say "正在读取最新高血压1256NOVO1000阴性样本检出率文件： $latest_control_detection_ratio";
    open CONTR,"$latest_control_detection_ratio" || die "can't open the $latest_control_detection_ratio!";
        %control_ditaction_ratio=map{chomp;(join"\t",(split/\t/,$_)[0..3]) => (join"\t",(split/\t/,$_)[4..6])}<CONTR>;
    close CONTR;
say "读取完成！！";
}

sub sub_case_ditacion_ratio{
    my $detection_ratio_prefix="snp_indel_detection_out";
    my $latest_detection_file=judge_lastest_sitedb2($case_detection,$detection_ratio_prefix);#查找最新的样本检出率文件
say '*'x25;
say "正在读取最新高血压检出率文件： $latest_detection_file";
    open CASEDI,"$latest_detection_file" || die "can't open the $latest_detection_file!";
        %case_ditacion_ratio=map{chomp;(join"\t",(split/\t/,$_)[0..3]) => (join"\t",(split/\t/,$_)[4..6])}<CASEDI>;
    close CASEDI;
say "读取完成！！";
}

########################## clinvar star ############################

sub read_clinvar_omim{
	my $file=shift;
	open(IN,$file)||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		$clivarPath{$t[0]}{path}=$t[6];
		$clivarPath{$t[0]}{star}=$t[10];
		$clivarPath{$t[0]}{OMIM}=$t[7];
	}
	close IN;
}

sub acmg_auto_judge{                    #wph add 190312       
	my $site_info = shift;
	my $sample    = shift;

	my @cut=split/\t/,$site_info;
	my @acmg=();
	my ($pep1,$pep2,$pep3)=("","","");
	my @site;

        if($cut[33] == 1 && $cut[8]=~/^nonsynonymous/ ){
                push @acmg,"PP3";
        }
        if($cut[34]!~/deleterious|damaging|-/i && $cut[8]=~/^nonsynonymous/ ){
                #push @acmg,"BP4";
        }

	my ($pep3_num1,$pep3_num2,$pep3_rate,$pep3_flag)=(0,0,0.0001,0);

        my $ref_alt=$cut[3]."_".$cut[4];
        my $repeat_site = $cut[15]."_".$cut[16];
        if( $cut[8]=~/^stoploss/ || ( $repeat_site eq "\._\." && $cut[8]=~/^nonframeshift/) ){
                #push @acmg,"PM4";
        }

	my $flag_BP7=0;
	if( $cut[8]=~/^synonymous/ || $cut[8]=~/^splicing/ )
	{
		#print $cut[28]."\n";
		$flag_BP7++ if($cut[28] ne "." && $cut[28]!~/SNV/ && $cut[28]!~/exon/ && $cut[28] > 0.6);
		$flag_BP7++ if($cut[29] ne "." && $cut[29]!~/SNV/ && $cut[29]!~/exon/ && $cut[29] > 0.6);
		$flag_BP7++ if($cut[32] ne "." && abs($cut[32]) >= 2 );
		if($flag_BP7==3 )
		{
			push @acmg,"PP3"; 
		}
	}

	my $flag_PS1=0;
	my $flag_truetype=0;
	my %star=('one'=>1,'two'=>2,'three'=>3,'four'=>4,'none'=>'0');
	my $clinvar_alle='';
	my ($truehgmd)=$cut[60]=~/genome:([^\|]+)/;
        my ($truecli)=$cut[61]=~/genome:([^\|]+)/;
        my ($truecliPep)=$cut[61]=~/pep:([^\|]+)/;
	if( $cut[8]=~/^nonsynonymous/ && $truehgmd!~/true_type/i  && $truecli!~/true_type/i && $truecliPep=~/true_type/)
	{
		@site=split /\|/,$cut[65];
		for(@site)
		{
			if( (split /:/,$_)[0] =~/true_type/i )
			{
				$pep2=(split /:/,$_)[3];
				$flag_PS1=1 if( $pep2=~/^Pathogenic/ || $pep2=~/^Likely pathogenic/ );
				$clinvar_alle= (split "\_", (split /:/,$_)[1] )[0];
				$flag_truetype=1 if( $star{ $clivarPath{$clinvar_alle}{star} } >=2 );
			}
		}
	}
	if($flag_PS1==1 && $flag_truetype==1 ){ push @acmg,"PS1"; }

	my $flag_PM5=0;
	$flag_truetype=0;
	my $flag_Ter=0;
	if( ($cut[8]=~/^nonsynonymous/) && $cut[63] ne "-")
	{
		@site=split /\|/,$cut[63];
		for(@site)
		{
			#next if( (split /:/,$_)[0] !~ /true_site/i );
			$pep1=(split /:/,$_)[4];
		#	$flag_PM5=1 if( $pep1 eq "DM" || $pep1 eq "DM?");
		}
	}
	if( ($cut[8]=~/^nonsynonymous/) && $cut[65] ne "-" ) 
	{
		@site=split /\|/,$cut[65];
		for(@site)
		{
			$flag_truetype=1 if( (split /:/,$_)[0] =~/true_type/i );
			next if( (split /:/,$_)[0] !~ /true_site/i );
			$pep2=(split /:/,$_)[3];
			#$flag_PM5=1 if( $pep2=~/^Pathogenic/ || $pep2=~/^Likely pathogenic/ );
			if( $pep2=~/^Pathogenic/ || $pep2=~/^Likely pathogenic/ )
                        {
                                $flag_PM5=1;
                                if( (split /:/,$_)[2]=~m/p\.(\S+)(\d+)(\S+)/ || (split /:/,$_)[2]=~m/p\.(\S)(\d+)(\*)/ )
                                {
                                        $flag_Ter=1 if($3=~"Ter" || $3 eq "*");
                                }
                        }
		}
	}
	if($flag_PM5==1 && $flag_truetype==0 && $flag_Ter==0 ){ push @acmg,"PM5"; }

	#my ($pep3_num1,$pep3_num2,$pep3_rate,$pep3_flag)=(0,0,0,0);
=head

	$pep3_flag=0;
	if( $cut[59] ne "-" )
	{
		if($cut[59] =~ /\((\d+)\/(\d+)\)/ )
		{
			$pep3_num1=$1;
			$pep3_num2=$2;
			$pep3_rate=sprintf("%.2f", $pep3_num1/$pep3_num2);
			$pep3_flag=1 if($pep3_num2>=7 && $pep3_rate>=0.6);
		}	
	}
	if($pep3_flag==1){push @acmg,"PM1";}

	my $flag_PVS=0;
	if($cut[8]=~/^splicing/)
	{
		@site=split /;/,$cut[7];
		for(@site)
		{
			$pep3=(split /:/,$_)[2];
			if( $pep3=~/c\.(\d+)([\+\-])(\d+)([ATCGN]+)>([ATCGN]+)/ )
			{
				if($3<3){ $flag_PVS = 1; }
			}
		}
		
	}
	if($cut[9] ne "." )
	{
		@site=split /,/,$cut[9];
		for(@site)
		{
			$pep3=(split /:/,$_)[3];
			if($pep3=~/c\.([ATCGN]+)(\d+)([ATCGN]+)/ )
			{
				if($2<4){ $flag_PVS = 1; }
			}
			elsif($pep3=~/c\.(\d+)del/ || $pep3=~/c\.(\d+)dup/)
			{
				if($1<4){ $flag_PVS = 1; }
			}
		}
	}
	#if($flag_PVS==1){ push @acmg,"PVS"; }

	$pep3_num1=0;
	$pep3_num2=0;
	$pep3_flag=0;

	if( $cut[8]=~/^stopgain/ || $cut[8]=~/^frameshift/ )
	{
		if( $cut[56] ne "-" )
		{
			if($cut[56] =~ /\((\d+)\/(\d+)\)/ )
			{
				$pep3_num1=$1;
				$pep3_num2=$2;
			#	$pep3_flag++ if($pep3_num1>0);
			}
			
		}
			
		if( $cut[57] ne "-" )
		{
			if($cut[57] =~ /\((\d+)\/(\d+)\)/ )
			{
				$pep3_num1=$1;
				$pep3_num2=$2;
				$pep3_flag++ if($pep3_num1>0);
			}
		}
		if( $pep3_flag>0 ){ $flag_PVS = 1; }
		else{
			push @acmg,"PM4";
		}
		#if($cut[8]=~/^frameshift/ && $pep3_flag==0 ){push @acmg,"PM4";}
	}

	if($flag_PVS==1){ push @acmg,"PVS"; }
=cut
	
	my $flag_PM2=0;
	my $db_dep_flag20=0;
	my @db_dep=split /\/\//, $cut[49];
	my ($db_dep_1, $db_dep_2, $db_dep_12, $db_dep_22);
	my ($db_dep_flag1,$db_dep_flag2) = (split /;/, $db_dep[2])[0,1];
	if($db_dep_flag1 eq "Y" )
	{
		$db_dep_1 = (split /\|/, $db_dep[0])[1];
		$db_dep_12= (split /,/, $db_dep_1)[0];
	}
	if($db_dep_flag2 eq "Y" )
	{
		$db_dep_2 = (split /\|/, $db_dep[1])[1];
		$db_dep_22= (split /,/, $db_dep_2)[0];
	}
	if($db_dep_flag1 eq "Y" && $db_dep_flag2 eq "Y" )
	{
		$db_dep_flag20=1 if( $db_dep_12>=20 && $db_dep_22>=20 );
	}
	if($db_dep_flag1 eq "Y" && $db_dep_flag2 eq "N")
	{
		$db_dep_flag20=1 if( $db_dep_12>=20);
	}
	if($db_dep_flag1 eq "N" && $db_dep_flag2 eq "Y")
	{
		$db_dep_flag20=1 if( $db_dep_22>=20);
	}

	if( ($cut[19] eq "." || $cut[19] eq "0") && ($cut[21] eq "." || $cut[21] eq "0") && $db_dep_flag20==1 )
	{
		$flag_PM2=1;
	}
	
	my $cont_ditacion_term=$control_ditaction_ratio{ join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-'); #NOVO1256样本检出率
	my $case_ditacion_term=$case_ditacion_ratio{     join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-');
	if( ( ( (split/\t/,$case_ditacion_term)[1]==1 && (split/\t/,$case_ditacion_term)[2] eq $sample ) || (split/\t/,$case_ditacion_term)[1]==0 )  && (split/\t/,$cont_ditacion_term)[1]==0 )
	{
		$flag_PM2=1 if( $cut[44] eq "snp");
	}
	if($flag_PM2==1){ push @acmg,"PM2"; }

	return @acmg;
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


