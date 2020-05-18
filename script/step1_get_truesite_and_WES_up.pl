use strict;
use warnings;
use Encode;
use utf8;
use autodie;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util qw/max min sum maxstr minstr shuffle/;
use feature qw/say/;
use lib '/PUBLIC/pipline/script/Report/XKFW/V3_180824/script';
use StandardClass;

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
	-cfg    样本信息文件
	-genelist 检测基因列表
	-samplelist 检测样本列表
END

my %disease;
my %database;
my %oldbase;
my %fubiaobasee;
my %database_2;
my %gene_disease_hash;
my %sub_sampleline;
my %clivarPath;
my %clivarPathsyn;
my %headOPT;
my %control_ditaction_ratio;
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
my $wesgene_info="/PUBLIC/pipline/database/siteFilterDB/product/WES/XK_IDT_subtype.info"; #wph add 190327
my $control_detetion="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/control_detection_ratio/170122_v0_1256wes_detection_ratio_out"; #zhouwei change 20170712 
my $clinvar_start_file="/PUBLIC/pipline/database/Clinvar/clivar_omim_stat.txt";
my $control_detection_path ||="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/control_detection_ratio";
my $case_detection ||="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/gxy_detection_ratio/";    #case集样本位点检出率


my $disease_script=judge_lastest_sitedb('Disease');           ###add panqi 190709
my %gene_disease_relation=&sub_omim_id_auto($disease_script); ###add panqi 190709
my $cusigene=$genelist;
my %gene=&genelist($cusigene);
my %detetion=&Deteion($control_detetion);

&sub_control_ditaction_ratio;
&sub_case_ditacion_ratio;
&read_clinvar_omim($clinvar_start_file);

my %samplelist=&readsamplelist($samplelist);
our %protein=(
	'A'=>'Ala','C'=>'Cys','D'=>'Asp','E'=>'Glu','F'=>'Phe','G'=>'Gly','H'=>'His','I'=>'Ile','K'=>'Lys','L'=>'Leu',
	'M'=>'Met','N'=>'Asn','P'=>'Pro','Q'=>'Gln','R'=>'Arg','S'=>'Ser','T'=>'Thr','V'=>'Val','W'=>'Trp','Y'=>'Tyr',
	'X'=>'*',
);  #wph add 190122

our @DatabaseHead=("ACMG对应条目","位点结论","疾病描述","基因描述","位点描述","参考文献","疾病名称","遗传模式","HGNC基因名","CDS","PEP","PDF_NMID");
      
##################################
open OUT, ">$out/$types.err.txt" || die "can't open the $out/$types.err.txt!";
open RT, "> $out/$types.right.txt" || die "can't open the $out/$types.right.txt!";
open RT2, "> $out/$types.right.2.txt" || die "can't open the $out/$types.right.2.txt!";#add right.2 file wph 18.8.21
if($type eq "IDT"){
	my $base=judge_lastest_sitedb($types);
	#my $fubase=judge_lastest_sitedb('IDTFB');
	%disease=read_disease($disease_file);
	%database=read_database($base,$type);
	#%fubiaobasee=read_database_2($fubase,$type);
	%fubiaobasee=read_database_2();
	%gene_disease_hash=read_genemiaoshu();
	&CS_get_basesite($path);
}else{
	print "$usage!"; 
}
close OUT;
close RT;
close RT2;


############################ cus filter site ###################################
sub CS_get_basesite{
	my $putin=shift;
	open IN, "$putin" || die "can't open the file!";
	my $count=1;
	my $headOUT="";
	my $headRT="";
	while(<IN>){
		my $prefix=(split/\./,(split/\//,$_)[-1])[0]; #file list
		chomp $_;
		my $file=$_;
		open FILE,"<:utf8","$file" or die "can't open the file!"; # file wangpenghui 181017
		chomp(my $head=<FILE>);
		if($count == 1){
			my ($namesHead,%head)=&FileHead($head);
			$headOUT= "#Flag\t$namesHead\n";
			$headRT="#Product\t$namesHead\t解读人\t末尾列\n";
			
			print OUT $headOUT;
			print RT $headRT;
		}
		while(<FILE>){
			chomp;
			my $disease_system='-';
			next if(/^CHR/i);
			my @tmp=split/\t/,$_;
			my $ID=join("\t",@tmp[0,1,3,4]);
			if(exists $detetion{$ID}){
				my ($ratio,$num)=split/\t/,$detetion{$ID};
				$tmp[73]=$num;#control 集检出人数
				$tmp[72]=$ratio; # control 集检出率
			}else{
				$tmp[73]="0"; #wph chaged 190103
				$tmp[72]="0";
			}
			if(exists $disease{$tmp[6]}){
				$disease_system=$disease{$tmp[6]};
			}else{
		                $disease_system="-";
			}
			my $re_info=join("\t",@tmp[0..68],$disease_system,@tmp[69..73]);
			my $Keeplabel=0;
			####
			my %currHash=();
			for(my $i=0;$i<@tmp;$i++){
				my $tit=$headOPT{$i};
				$currHash{$tit}=$tmp[$i];
			}
			################
			my($flag_exists,$flag_DP,$flag_Rate,$flag_protein,$flag_dm,$flag_benign,$flag_hot,$flag_Negative,$mu_type,$flag_snp,$flag_Pathsyn,$errorTag)=StandardClass::Standard_confirm(\%currHash,\%gene,$prefix);
			     #基因在基因列表，深度满足，频率满足，蛋白预测满足，至少一个有害，全良性，在热点区，删除点， 其余点1 exon等0 异常标签
			my ($databaseDP,$acmg)=StandardClass::acmg_auto_judge(\%currHash,$prefix,$re_info,\%clivarPath,$flag_Rate,\%control_ditaction_ratio,\%case_ditacion_ratio);
			$errorTag.=";BaseDP<20" if($databaseDP==0);
			my $pick_num=@$acmg;
			my $acmgLine=join(";",@$acmg);
			$errorTag="$errorTag;$acmgLine"if ($pick_num>0);
			if( $flag_DP==1 && $flag_Negative==0) #前提：在分析列表，深度满足，非待删除点
			{
				if($flag_hot==1){
					$Keeplabel=1;  #热点区无条件满足
					$errorTag.=";SVI-rev";
				}else{
					if(($flag_Rate==1 && $databaseDP==1)  &&  ($flag_protein==1))
					{       
						$Keeplabel=1; #满足频率和蛋白预测，则保留
						if($currHash{"ExonicFunc"} =~/^synonymous/){
							if($flag_dm==0){
								$Keeplabel=0;
								if($flag_Pathsyn==1 && ($acmgLine=~/PM2/ || $acmgLine=~/PP3/)){
									$Keeplabel=1;
								}
							}
						}elsif($pick_num==0 && $flag_snp==1 ){ #如果没有被捡回的条目，满足下面则删除,indel不捡回
							$Keeplabel=0 if($mu_type==0 && $flag_benign==1 || $mu_type==1 && $flag_dm==0 );
						}
					}else{
						if ($pick_num>0 && $flag_snp==1 ){  
							if($currHash{"ExonicFunc"} !~/^synonymous/){
								$Keeplabel=1;	
							}else{
								$Keeplabel=1 if($flag_Pathsyn==1 && ($acmgLine=~/PM2/ || $acmgLine=~/PP3/));
							}
						}
						if($flag_snp==0){
						 	$Keeplabel=1 if($errorTag!~/AF/ );
						}
					}
				}
			}
			&right_site($prefix,$re_info,$flag_exists,$Keeplabel,$errorTag); 
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

###########sub file head###
sub FileHead{
	my $head=shift;
	my $head_old=$head;
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
	my $headnew="ID\t$names\tnone\t$database{head}\tcase_ratio\tcase_num\tcase_ID";
	
	for (my $i=0;$i<@head_name;$i++) {
		$headOPT{$i}=$head_name[$i];
	}
	return $headnew;
}


############
sub right_site{
	my @arr=@_;
	my $prefix=$arr[0];
	my $all=$arr[1];
	my $falg_exsits=$arr[2];
	my $flag_RT=$arr[3];
	my $flag_wrong=$arr[4];
	
	($all,my @tmp3)=&add_sitebase_info($all);
	my @tmp=split/\t/,$all;
	my $key=join("_",@tmp[0,1,3,4]);
	my $key22=join("\t",@tmp[0,1,3,4]);

	$RTsampleIDs{$prefix}=1 if($flag_RT==1 && $falg_exsits==1);
	$dupsites{$key22}{$prefix}=1 if($flag_RT==1 && $falg_exsits==1);
	my @allg=split /;/,$tmp[6];

	if(!exists $case_ditacion_ratio{$key22}){$case_ditacion_ratio{$key22}="0\t0\t-";} #wph chaged 190103
	my %outhash;
	my @resultLine=();
	my %sitetime=();
	if(exists $database{$key} ){  ##zhouwei  如果是false位点的话，那么解读库中的HGMD id号要写成-，否则不会进行匹配。在right文件中就会被过滤掉
		for(@{$database{$key}}){
			my @ttt=split /\t/,$_;
                        my $key2="$key\_$ttt[6]";
			if(!exists $sitetime{$key2}  || $sitetime{$key2} < $ttt[-1]){
				$sitetime{$key2}=$ttt[-1];
			}
		}
		for(@{$database{$key}}){
			my @ttt=split /\t/,$_;
			my $key2="$key\_$ttt[6]";
			$outhash{$key2}=1; 
			my $currentLind="";
        		my $currentTag="";
			my $databaseline="";
			foreach my $k(@allg){
				if($tmp[70]=~/^chr/ && exists $gene_disease_relation{$k}{$ttt[6]} ){  ###add panqi 190709
					$tmp[70]=$gene_disease_relation{$k}{$ttt[6]};
				}
			}
			my $all2=join("\t",@tmp);
			if($ttt[11] !~ /NM/){
				$ttt[11] = $tmp[12];
			}
			if($sitetime{$key2} eq  $ttt[-1]){
				$databaseline=join("\t",'-',@ttt[0..11]);
				$currentTag="Done$ttt[12]";
				$currentLind="$all2\t$databaseline";
				push @resultLine,[$currentLind,$currentTag];
			}
		}
	}

	if(exists $fubiaobasee{$key}){
		for(@{$fubiaobasee{$key}}){
			my $line2=$_;
			my @aaa=split /\t/,$line2;
			my $key2="$key\_$aaa[6]"; 
			my $currentLind="";
                        my $currentTag="";
                        my $databaseline="";
			foreach my $k(@allg){
				if($tmp[70]=~/^chr/ && exists $gene_disease_relation{$k}{$aaa[6]} ){  ###add panqi 190709
					$tmp[70]=$gene_disease_relation{$k}{$aaa[6]};
				}
			}
			my $all2=join("\t",@tmp);

			if(!exists $outhash{$key2}){
				if($aaa[11] !~ /NM/){
					$aaa[11] = $tmp[12];
				}
				$databaseline=join("\t",'-',@aaa[0..11]);
				$currentTag="Done$aaa[12]";
				$currentLind="$all2\t$databaseline";
				push @resultLine,[$currentLind,$currentTag]; 
				$outhash{$key2}=1;
			}
		}
	}
			

	for(@tmp3){
		my @tmp2=split /\t/,$_;
		my $databaseline="-\t".$_."\t$tmp[12]";
		my $key2="$key\_$tmp2[6]"; 
		if(!exists $outhash{$key2} ){#wph changed
			foreach my $k(@allg){
				if($tmp[70]=~/^chr/ && exists $gene_disease_relation{$k}{$tmp2[6]} ){  ###add panqi 190709
					$tmp[70]=$gene_disease_relation{$k}{$tmp2[6]};
				}				
			}
			my $all2=join("\t",@tmp);
			my $currentLind="";
                        my $currentTag="";
			if(!exists $outhash{$key2}){
				$currentTag="Todo";
				$currentLind="$all2\t$databaseline";
				push @resultLine,[$currentLind,$currentTag]; 
				$outhash{$key2}=1;
			}

		}
	}

	$flag_wrong=";$flag_wrong" if($flag_wrong!~/^;/);
	foreach my $p (@resultLine) {
		print RT "$p->[1]$flag_wrong\t$prefix\t$p->[0]\t$case_ditacion_ratio{$key22}\n"  if($flag_RT==1 && $falg_exsits==1);
		print OUT "$p->[1]$flag_wrong\t$prefix\t$p->[0]\t$case_ditacion_ratio{$key22}\n" if($flag_RT==0 && $falg_exsits==1);
		print RT2 "$p->[1]$flag_wrong\t$prefix\t$p->[0]\t$case_ditacion_ratio{$key22}\n" ;
	}
}

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
##########读数据库
sub read_disease{
	my $disease=shift;
	my %h;
	die "$disease is not exists" unless -e $disease;
	open IN,"<:utf8", $disease or die "Can not open file $disease";
	while(<IN>){
		chomp;
		my ($gene,$system)=split/\t/,$_,2;
		$h{$gene}=$system;
	}
	close IN;
	return %h;
}
sub read_genemiaoshu{#wph add 190327
	my %gene_disease_hash=();
	open INN,"<:utf8",  $wesgene_info or die $!; #"<:utf8"
	while(<INN>)
	{
		chomp;
		next if($_=~/疾病英文/);
		my @cut=split /\t/,$_;
		push @{$gene_disease_hash{$cut[1]}}, join("\t",@cut[7,6,3,4,1]);
	}
	close INN;
	return %gene_disease_hash;
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
		if($t[2]=~/p\./){
			my ($gene,$Pref,$Palt)=$t[2]=~/\((.*)\):.*\(p\.([A-Za-z]+)\d+([A-Za-z=]+)\)/; #NM_000410.3(HFE):c.193A>T (p.Ser65Cys)
			if($Palt && $Pref && $Palt eq $Pref || $Palt && $Palt eq '='){
				$clivarPathsyn{$gene}=1 if($t[6] =~/Pathogenic/ || $t[6] =~/Likely pathogenic/);
			}
		}
	}
	close IN;
}

################################## site database #############################################
sub read_database{
	my $shuju=shift;
	my $product=shift;
	open IN, "<:utf8", "$shuju" || die "can't open the file";
	my %hash;
	$hash{'head'}=(join"\t",@DatabaseHead); 
	my %heads=();
	while(<IN>){
		chomp;
		my @tmp=split/\t/,$_;
		if(/^#flag/){
			for(my $i=0;$i<@tmp;$i++) {
				$heads{$tmp[$i]}=$i; 
			}
		}else{
			my @keys=();
			my @keyarry=("CHROM","Start","REF","ALT");
			for(my $i=0;$i<@keyarry;$i++){
				my $index=$heads{$keyarry[$i]};
				if($keyarry[$i] eq 'REF' && $tmp[$index] eq '.'){
					$tmp[$index]='-';
				}
				if($keyarry[$i]  eq 'ALT' && $tmp[$index] eq '.'){
					$tmp[$index]='-';
				} 
				push @keys,$tmp[$index];
			}
			my $key=join("_",@keys); 
			my $disnameindex=$heads{'疾病名称'};
			my $key2="$key\_$tmp[$disnameindex]"; 
			my $n=@tmp;
			my $NMindex=$heads{'PDF_NMID'};
			if($tmp[$NMindex] !~/NM/){
				$tmp[$NMindex]='-';
			}
			if($tmp[0]=~/^T1/){
				$sub_sampleline{$key2}=1;
			}
			my @a=();
			for(my $i=0;$i<@DatabaseHead;$i++) {
				my $index=$heads{$DatabaseHead[$i]};
				push @a,$tmp[$index];
			}
			push @a,$tmp[0];
			push @a,$tmp[-1];
			my $info = join("\t",@a);
			push @{$hash{$key}},$info ;
		}	
	}
	close IN;
	return %hash;
}

sub read_database_2{ ################### sunxq 20180408
	#my $shuju=shift;
	#my $product=shift;
	my $shuju="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_zhufu/WES_zhubiao_fubiao.txt";
	open IN, "<:utf8", "$shuju" || die "can't open the file";
	my %hash;
	my $aa=<IN>;
	chomp $aa;
	$hash{'head'}=join"\t",(split /\t/,$aa)[79..90];
	while(<IN>){
		chomp;
		my @tmp=split/\t/,$_;
		my $key=join("_",@tmp[3,4,6,7]);
		$tmp[90]='-' if(!exists $tmp[90]);
		my $info=join("\t",@tmp[79..90],$tmp[0]);
		push @{$hash{$key}},$info;
	}
	close IN;
	return %hash;
}




####################
sub genelist{
	my $file=shift;
	my %hash;
	my %sample_panel;	
	my $sample;
	my @genes;#wangph change 20180606
	my $ge;#wph
	my $disease;
	my @empty;
	my @gxy_disease=split /\n/, `cat /PUBLIC/pipline/database/siteFilterDB/product/GXY/单基因高血压疾病`;
	my @cus_disease=split /\n/, `cat  /PUBLIC/pipline/database/siteFilterDB/product/CuS/心源性猝死`;
	my @dm_disease=split /\n/, `cat /PUBLIC/pipline/database/siteFilterDB/product/DM/单基因糖尿病`;
	print "gene $file\n";
	open IN,"<:utf8","$file" || die "can't open the file!";
	while(<IN>){
		chomp;
		$sample=(split/\t/,$_)[0];#wph changed 20180615
		$disease=(split/\t/,$_)[1];#wph changed 20181207
		if($disease=~/单基因高血压|心源性猝死|单基因糖尿病/){#wph changed 20190327
			$sample_panel{$sample}=
			($disease eq "单基因高血压")		? @gxy_disease	:
			($disease eq "心源性猝死")		? @cus_disease	:
			($disease eq "单基因糖尿病")		? @dm_disease	:  @empty;
			@genes=split/, /,(split/\t/,$_)[2];
			foreach $ge (@genes){
				$hash{$sample}{$ge}=0;
			}
		}else{
			@genes=split/, /,(split/\t/,$_)[2];#wph 180606
			my $match = grep { $_ =~ /$disease/ } $sample_panel{$sample};
        		foreach $ge (@genes){
				if(not $match and !exists $hash{$sample}{$ge}){
					$hash{$sample}{$ge}=1;
				}elsif(exists $hash{$sample}{$ge}) {
					next;
				}else{
            				$hash{$sample}{$ge}=1;
				}
			}
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

sub sub_omim_id_auto{                                        ###add panqi 190709
	my $disease_script=shift;
	my %gene_disease_relation=();
	my @cut;
	my @arr;
	my $disease_name;
	my $tmp2;

	open INN,"<:utf8", $disease_script or die $!;
	while(<INN>){
		chomp;
		next if($_=~/^#NewID/);
		@cut=split /\t/,$_;
		$disease_name= (split /\【/,$cut[4])[0]; 		
		$gene_disease_relation{$cut[3]}{$disease_name}=$cut[0];
		@arr=split /\|/,$cut[12] if($cut[12] ne "-" && $cut[12] ne "NULL");
		next if($cut[12] eq "-" || $cut[12] eq "NULL");
		for(@arr)
		{
			$disease_name=(split /\_/, $_)[1];
			$gene_disease_relation{$cut[3]}{$disease_name}=$cut[0];
		}
	}
	close INN;
	return %gene_disease_relation;
}

#############自动添加注释信息
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
	if( ($cut[5]=~/splicing/i and $cut[5]!~/exon/) or $cut[8]=~/splicing/ or $cut[8]=~/^synonymous SNV/i ){
		$splicingds="该变异经可变剪切软件进行预测，";
        if( $cut[28] ne "." && $cut[29] ne "." && $cut[32] ne "." && $cut[28] > 0.6 && $cut[29] > 0.6 && abs($cut[32]) >=2)
        {
            $splicingds.="dbscSNV1.1、SPIDEX软件预测该变异对基因或基因产物有害。";
        }
        elsif( $cut[28] ne "." and $cut[29] ne "." and $cut[29] !~/ENSG/ and  $cut[29] !~/SNV/ and $cut[28] > 0.6 and $cut[29] > 0.6)
        {
            $splicingds.="dbscSNV1.1软件预测该变异对基因或基因产物有害。";
        }
        elsif($cut[32] ne "." and abs($cut[32]) >=2)
        {
            $splicingds.="SPIDEX软件预测该变异对基因或基因产物有害。";
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

	my $software_pred='';
	my @softwares=split /;/,$cut[34];
	my @ssss;
	my $sss_n=0;
	if($softwares[0] =~/Deleter/i || $softwares[0] =~/dama/i){push @ssss,'SIFT';$sss_n++;}  #SIFT、Polyphen2_HVAR、Polyphen2_HDIV、M-CAP
	if($softwares[1] =~/Deleter/i || $softwares[1] =~/dama/i){push @ssss,'Polyphen2_HVAR';$sss_n++;}
	if($softwares[2] =~/Deleter/i || $softwares[2] =~/dama/i){push @ssss,'Polyphen2_HDIV';$sss_n++;}
	if($softwares[3] =~/Deleter/i || $softwares[3] =~/dama/i){push @ssss,'M-CAP';$sss_n++;}
	my $sssware=join("、",@ssss);
	if($sss_n>0){
		$software_pred=$sssware."软件预测该变异对基因或基因产物有害。"
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
		my @tmp=split /\t/, $jiyinmiaoshu[$i];
		@add_info[0,1]=($cut[43]=~/true/i) ? ("待解读","待解读") : ("自动","临床意义未明") ; #ACMG条目,#致病性结论
		@add_info[2,3,6,7,8]= @tmp[0,1,2,3,4];  #致病基因 疾病名称 疾病描述 遗传模式 临床指导建议 基因描述    临床表现 疾病大类
		@add_info[9,10]=@cut[13,14]; 
		($add_info[10],$extra_type)=&pep_trans($cut[14]); #CDS    PEP
		#if($mutation_type ne $extra_type){
		#	print "$mutation_type:$cut[8]_$cut[13]_$cut[14]"."|".$extra_type."_$add_info[11]\n";
		#}	
		if($mutation_type eq "?" or $cut[8] =~/unknown/i){  #wph add 190418)
			$mutation_type=$extra_type;
			$cut[8] =
			($extra_type eq "错义")					? "nonsynonymous SNV"	    :
			($extra_type eq "同义")					? "synonymous SNV"	    :
			($extra_type eq "无义")					? "stopgain"		    :
			($extra_type eq "终止缺失")				? "stoploss"		    :
($extra_type eq "剪切" && $cut[13] !~/del|->|>-/i && $cut[5] ne "intronic")	? "splicing SNV"	    :
			($extra_type eq "剪切" && $cut[13] =~/(del|->|>-)/i)	? "splicing INDEL"	    :
			($extra_type eq "非移码" && $cut[13] =~/del/) 		? "nonframeshift deletion"  :
			($extra_type eq "非移码" && $cut[13] =~/(dup|ins)/) 	? "nonframeshift insertion" :
			($extra_type eq "移码" && $cut[13] =~/del/)		? "frameshift deletion"	    :
			($extra_type eq "移码" && $cut[13] =~/(dup|ins)/)	? "frameshift insertion"    : "SNV";
		}
		my $muty="";#wph changed 1910104 $cut[14]${mutation_type}";
		if ($cut[5]=~/splicing/ and $cut[5]!~/exon/ or $extra_type eq "剪切"){
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
		$muty=~s/^-//g;
		$add_info[4]="本筛查检测出$add_info[8]基因的变异，变异位点为$cut[13]，查询ClinVar等公共数据库显示，$cut[13]变异${muty}。$ratio_detection"."${software_pred}${splicingds}${wenxianzhichi}行ACMG标准，判定该变异的致病性为临床意义未明。具体情况请结合临床相关资料综合判断。";#位点描述
		#$add_info[6]="-"; #参考文献
		push @add_info_all, join("\t", @add_info[0..10]);
	}
	#$add_info_all=join("\t", @add_info[0..11]);
	$site_info=join("\t",@cut);
	my $len=@cut;
	return $site_info,@add_info_all;
}

sub get_jiyinmiaoshu{ #wph add 190327
	my $gene=shift;
	my $cha_gene=$gene=~s/\(.*\)//gr;
	my @all_gene=split/;/,$cha_gene;
	my %hash=();
	my @cut;
	my @all_term;

	foreach my $eachgene (@all_gene){
		if(exists $gene_disease_hash{$eachgene}){
			foreach my $item ( @{$gene_disease_hash{$eachgene}} ){
				push @all_term, $item;  
			}
		}
		else{
			push @all_term,"-\t-\t-\t-\t$eachgene";
		}
	}
	return @all_term;
}

sub pep_trans{ #wph add 190122 190402 add type
	my $p=shift;
	my $type="";
	if($p=~/p\.([A-Y]+)(\d+)([A-Y]+)(.*)/g){
		my $locaa=$2;
		my $firstaa =join"",(map{$protein{$_}//'-'}(split//,$1));
		my $secondaa=join"",(map{$protein{$_}//'-'}(split//,$3));
		my $lastaa=$4;
		$p="p\.${firstaa}${locaa}${secondaa}${lastaa}";
		if($p =~/fs/){
			$type="移码";
		}elsif($firstaa eq $secondaa){
			$type="同义";
		}elsif($secondaa =~ /\*/){
			$type="无义";
		}elsif($firstaa =~ /\*/){
			$type="终止缺失";
		}else{
			$type="错义";
		}
	}elsif($p=~/p\.([A-Y]+)(\d+)fs/g){
		my $locaa=$2;
		my $firstaa =join"",(map{$protein{$_}//'-'}(split//,$1));
		$p="p\.${firstaa}${locaa}fs";
		$type="移码";
	}elsif($p eq "-"){
		$type="剪切";
	}elsif($p =~/fs/i){
		$type="移码";
	}elsif($p =~/X/i){
		$type="无义";
	}else{  
		$type="非移码";
		return $p,$type;
	}
	return $p,$type;
}

#################判断最新数据库地址
sub judge_lastest_sitedb{
	my $product_type=shift;
	my %hash_product_sitepath=();
	my @all_date;
	my $lastest_date;
	my $lastest_sitedb_path;
	my $gxy_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY"; #zhouweichange 20170712
	my $cusi_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_CuS"; #zhouweichange 20170629
	my $wes_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES";
	my $wes_sampledb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_report_base/";
	my $wes_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_fubiao/";
	my $genetic_sitedb_path="/ifs/TJPROJ3/HEALTH/zhongying/00research/00report/01database";
	my $cr_sitedb_path="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/CR";  #zhouweichange 20170629
	my $NNY_sitedb_path="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/Nuonanyu";  #zhouweichange 20170629
	my $genetic_black_box="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/genetic_heihe";   #zhouweichange 20170629
	my $cusi_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Fubiao_CuS/";  #sunxiaoqing
	my $gene_descript_path="/PUBLIC/pipline/database/Auto_pipline/gene_description/";#wangpenghui        ###add panqi 190709
	my $disease_descript_path="/PUBLIC/pipline/database/Auto_pipline/disease_description/";#wangpenghui  ###add panqi 190709
	$hash_product_sitepath{'FWGXY'}=$gxy_sitedb_path;
	$hash_product_sitepath{'GXY'}=$gxy_sitedb_path;
	$hash_product_sitepath{'CS'}=$cusi_sitedb_path;
	$hash_product_sitepath{'IDT'}=$wes_sampledb_path;
	$hash_product_sitepath{'IDTFB'}=$wes_sitedb_path_fubiao;
	$hash_product_sitepath{'YC'}=$genetic_sitedb_path;
	$hash_product_sitepath{'YC133'}=$genetic_sitedb_path;
	$hash_product_sitepath{'CR'}=$cr_sitedb_path;
	$hash_product_sitepath{'NNY'}=$NNY_sitedb_path;
	$hash_product_sitepath{'black'}=$genetic_black_box;
	$hash_product_sitepath{'fubiao'}=$cusi_sitedb_path_fubiao;
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
