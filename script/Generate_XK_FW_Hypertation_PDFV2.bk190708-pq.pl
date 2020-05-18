#!/usr/bin/perl
#use strict;
use warnings;
use feature qw/say/;
use File::Find;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw/max min sum maxstr minstr shuffle/;
use Getopt::Long;

my ($final_2_path,$step,$out,$public,$case_detection,$hyper_sample_info_path,$gene_desease_path,$control_detection_path,$disease_file,$disease_ratio_file,$latest_gxybase_path,$latest_gxybase_path_fubiao,$latest_gxyreportbase_path,$latest_deletion_path,$hypers_cfg_addsex_path,$NOVO_1256_AF,$HYper_new_AF,$latest_gxybase_old_path);
GetOptions(
    'final_2_list=s' => \$final_2_path, #必选	
    'step=s' => \$step,  ##可选  （1,2,或者不选）
    'out=s' => \$out, #必选
    'public=s'=> \$public, #可选
    'case_detection=s'=> \$case_detection, #可选
    'control_detection=s'=> \$control_detection_path, #可选
    'doctor_patient=s'=> \$hyper_sample_info_path, #可选
    'gene_disease_discpt=s'=> \$gene_desease_path, #可选
    'gene_disease_system=s'=> \$disease_file, #可选
    'latest_gxybase_path=s'=> \$latest_gxybase_path, #可选
    'latest_gxyreportbase_path=s' =>  \$latest_gxyreportbase_path,#可选
    'latest_gxybase_old_path=s' => \$latest_gxybase_old_path, #可选
    'latest_deletion_path=s' => \$latest_deletion_path, #可选
    'add_sex_config=s'=> \$hypers_cfg_addsex_path, #可选
    'novo_1256_af=s'=> \$NOVO_1256_AF, #可选
    'hype_last_af=s'=> \$HYper_new_AF, #可选

);

###########################默认文件路径设置######################################
$public ||="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/"; #高血压公共数据库路径
$case_detection ||="$public/gxy_detection_ratio";    #高血压case集样本位点检出率
$control_detection_path ||="$public/control_detection_ratio";  #高血压control集样本位点检出率  
$hyper_sample_info_path ||="$public/gxy_sample_info";  #医生患者信息单
$gene_desease_path ||="$public/gene_deseases";   #高血压基因疾病描述
$disease_file ||="/PUBLIC/pipline/database/data/gene/gene_add_desase.xls";  #高血压基因及对应的疾病系统
$disease_ratio_file ||="/PUBLIC/pipline/script/Report/XKFW/V3_180824/script/GXY_incidence.txt";
$latest_gxybase_path ||="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY";  #高血压位点库 #zhouweichange 20170629
$latest_gxybase_path_fubiao ||= "/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY_fubiao/"; 
$latest_gxyreportbase_path ||="/PUBLIC/pipline/database/knowledge_database/sitedatabase/GXY_report_base/"; #高血压出报告位点库 #sunxiaoqing 180511
$latest_gxybase_old_path ||= "/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY//201711221646/XK_gxybase201711221646.right.txt";  #高血压位点库 #sunxqchange 20171130
$latest_deletion_path ||= "/PUBLIC/pipline/database/knowledge_database/sitedatabase/filter_site"; # 待删除库
$hypers_cfg_addsex_path ||="$public/gxy_all_cfg_addsex"; #实时更新的高血压config文件，和对应的gxy_detection_ratio对应的人数一致，如果不一致，请检查mergvcf是否执行成功
$NOVO_1256_AF ||="$public/AF/1256WES_Novogene2016_AF.txt"; #novo control 1000WES等位基因频率
$HYper_new_AF ||="$public/AF/final_head_site_all.frq.out.base.right"; #高血压累积case等位基因频率，目前未实时更新，需要今后完善
my $clinvar_start_file="/PUBLIC/work/sunxiaoqing/anno_updata/database/Pathogenicity_Check/clivar_omim_stat.txt";

#############################主脚本区域###########################################
my (%case_ditacion_ratio, %control_ditaction_ratio, %doctor_patient_info, %genelist, %gene_disease, %disease, %disease_ratio, %hash_base, %hash_base2 ,%hash_base_fubiao, %hash_base_fubiao2 ,%database_site ,%hash_base_2,%hash_del,%right_old, %sampleid_sex, %novo_1256_af, %hyper_last_af, %clivarPath);

our %protein=(
	'A'=>'Ala','C'=>'Cys','D'=>'Asp','E'=>'Glu','F'=>'Phe','G'=>'Gly','H'=>'His','I'=>'Ile','K'=>'Lys','L'=>'Leu',
	'M'=>'Met','N'=>'Asn','P'=>'Pro','Q'=>'Gln','R'=>'Arg','S'=>'Ser','T'=>'Thr','V'=>'Val','W'=>'Trp','Y'=>'Tyr',
	'X'=>'*',
);  #wph add 190122

my @titel=(
"Sample",
"CHROM\tStart\tEnd\tREF\tALT\tFunc\tGene\tGeneDetail\tExonicFunc\tAAChange\tNo_annovar_trans\ttransvar_Aachange\tpdf_NM_ID\tpdf_CDS\tpdf_PEP\tgenomicSuperDups\trepeat\tavsnp150\tRare_SNP\tgnomAD_exome_EAS\t1000g2015aug_eas\texac2015_03newEAS\tCombineWES\tphyloP20way_mammalian\tphyloP20way_mammalian_rankscore\tphastCons20way_mammalian\tphastCons20way_mammalian_rankscore\tInterpro_domain\tada_score\trf_score\tHET/HOM\tdpsi_max_tissue\tdpsi_zscore\tProtein\tSIFT;Polyphen2_HVAR;Polyphen2_HDIV;M-CAP;Revel\tSIFT\tPolyphen2_HDIV\tPolyphen2_HVAR\tMCAP\tMutationTaster/CADD_phred\tcytoBand\tGENE_strand\tWarning\tTYPE\tsnp/indel\tDP\tALT_rate\tcds_stop_dis\tbed_edge_dis\tgnomad_depth\tuniprot_domain\tgene_disease\tsite_hotspot\tOMIM\tgene_HGMD_pathogen_contex\tgene_Clivar_pathogen_contex\tHGMD_pathogen_PVS1\tClivar_pathogen_PVS1\tHGMD_pathogen_15bp\tClivar_pathogen_15bp\tTrue_Flase_HGMD\tTrue_Flase_Clinvar\ttrue_genome_hgmd_contex\ttrue_pep_hgmd_contex\ttrue_genome_clivar_contex\ttrue_pep_clivar_contex\tVCF_INFO\tVCF_FORMAT\tVCF_sampleID\tchr_hg38\ts_hg38\te_hg38\tref_hg38\talt_hg38",
#"CHROM\tStart\tEnd\tREF\tALT\tFunc\tGene\tGeneDetail\tExonicFunc\tCDS\tPEP\tAAChange\tcytoBand\tavsnp150\tRare_SNP\tgnomAD_genome_ALL\tgnomAD_exome_ALL\tgnomAD_exome_EAS\tgnomAD_exome_SAS\t1000g2015aug_all\t1000g2015_Chinese\t1000g2015_exac03newALL\t1000g2015_exac03newEAS\t1000g2015_CombineWGS\t1000g2016_CombineWES\tDP\tDP4\tDP4/DP\tALT_RATE\tMQ\tHET/HOM\tWarning\tStrand prefrence\tProtein\tSIFT;Polyphen2_HVAR;Polyphen2_HDIV;M-CAP;Revel\tInterpro_domain\tphyloP20way_mammalian\tphyloP20way_mammalian_rankscore\tphastCons20way_mammalian\tphastCons20way_mammalian_rankscore\tada_score\trf_score\tdpsi_zscore\tTYPE\tsnp/indel\tHGMD;COSMIC;Clinvar\tAA_type\tcosmic81\tclinvar_type\tchr_DB\ts_DB\te_DB\tGene_DB\tHGMD_ID_DB\tAachange_DB\tDisease_name\tclass\tChina/English/inheritance\tpaper\tPubMed\tPhenotypes\tpubmed_id\tMim_Number\tGene_Symbols\tomim_id\tsymbol\tname\talias_symbol\tccds_id\tchr_hg38\ts_hg38\te_hg38\tref_hg38\talt_hg38",
"Disease_system",
"ACMG条目\t致病性结论\t疾病描述\t基因描述\t位点描述\t参考文献\t疾病名称\t遗传模式\t基因名称\tCDS\tPEP\t临床指导建议\t临床表现\tPDF_NMID",
"疾病大类\t疾病亚类",
"Sample_name\tDotcor\tHospital\tSample_gender",
"All_Sample_name\tAll_Dotcor\tAll_Hospital\tAll_Sample_gender",
"control_ditaction_ratio\tcontrol_ditaction_number\tcontrol_ditaction_sample_ID\tcontrol_ditaction_all_number",
"case_ditaction_ratio\tcase_ditaction_number\tcase_ditaction_sample_ID\tcase_ditaction_all_number",
"case_control_rates",
"control_N_ALLELES\tcontrol_SAMPLE_NUMBER\tcontrol_REF_AF\tcontrol_ALT_AF",
"case_N_ALLELES\tcase_SAMPLE_NUMBER\tcase_REF_AF\tcase_ALT_AF\t解读人\t末尾列",
);

#*************#
#执行代码*****#
&alljob;      #
#END**********#
#*************#

sub alljob{
say '*'x25;
say "诺禾心康和阜外高血压模块：";
say "step参数可选1和2，不选或者其他参数均为执行所有程序！";
say '*'x25;

    (!$step)   ? (&sub_step_1 and &sub_step_2) :
    ($step==1) ? &sub_step_1 :
    ($step==2) ? &sub_step_2 : (&sub_step_1 and &sub_step_2);
}

sub sub_step_1{
	&sub_doctor_patient_info;
	&sub_gene;
	&sub_gene_disease;
	&sub_disease;
	&sub_disease_ratio;
	&sub_hash_base;
	&sub_hash_fubbase;
	&sub_hash_oldbase;
	&sub_sampleid_sex;
	&sub_case_ditacion_ratio;
	&sub_control_ditaction_ratio;
	&sub_novo_1256_af;
	#&sub_hyper_last_af;
	&read_clinvar_omim($clinvar_start_file);


	open ALL,  ">$out/GXY_sample_all_info.txt"   or die $!;
	open RIGHT,">$out/GXY_sample_rigth_info.txt" or die $!;
	open WRONG,">$out/GXY_sample_wrong_info.txt" or die $!;
#open READ ,">$out/GXY_sample_site_database_read.txt" or die $!; ##changed sunxq 20171121
say WRONG <<'END';
高血压位点过滤Flag类型：
DP<4：Depth深度小于4；
AF1000_4>0.01：4列千人频率值有存在大于0.01；
No_disease：表示这个位点的基因未找到对应的疾病；
Man_pregnancy_hypertension：男性过滤掉妊娠高血压；
Number_case>100：case集的位点人数大于100均被过滤掉；
END
#DP4_ratio<0.1：DF4/DP小于0.1；
#Number_control>3：control集的位点人数大于3均被过滤掉；
#protoin_predict_filter：SIFT,Polyphen,M-cap预测没有匹配到“deleterious|damaging|-”，并其中不包括“.|-”;
#true_positive位点：不做过滤，前提满足Man_pregnancy_hypertension;

	say ALL   join"\t",@titel;
	say RIGHT join"\t",("Flag",@titel);
	say WRONG join"\t",("Flag",@titel);
	#say READ  join"\t",("-\t-\t-",@titel[0..3]);
	say "注意：对于所有位点的key均采用 chrom start ref alt，而一些数据库的INDEL位点的ref和alt表示方法不同，从而导致INDEL位点出现部分数据为----";
	my %sample_number;
	open FI,"$final_2_path" || die "can't open the $final_2_path!";
	while(<FI>){
		chomp;
		next if(/^$/);
		(say "此文件不存在！请检查文件路径是否正确！！\n$_" and next) unless (-e $_);
		(say "此文件存在，但是为空！请检查该文件生成是否正确！！\n$_" and next) if (-z $_);

		my $name=basename($_);
		my $sample_id=(split/\./,$name)[0];
		++$sample_number{$sample_id};
		open IN,$_ || die "can't open the $_!";
		while(<IN>){
			next if(/^CHROM/ || /^Chr/);
			chomp;
			my @cut=split/\t/;
			my $key=join"\t",@cut[0,1,3,4];
			my $cha_gene=$cut[6]=~s/\(.*\)//gr; 
			my @all_gene=split/;/,$cha_gene;  ###changed by sunxq 20171121
			##################
			my $key2;
			if($cut[43]=~/true/i){
					  $key2=join("\t",@cut[0,1,3,4]); #add end HGMD id 2016.08.23
			}else{
				if($cut[44]=~/snp/i){     # add by sxq 2017.8.15
						$key2=join("\t",@cut[0,1,3,4]);
				}else{
						$key2=join("\t",@cut[0,1,3,4]);
				}
			}

			if($cut[0] eq "16" and $cut[1]==1258036 )
			{
				$cut[1]=1258036;
			}

			if(exists $database_site{$key2} && $database_site{$key2}!~/UNKNOWN/i   && ($cut[9] =~/UNKNOWN/ || $cut[9] eq './/.,')){
				my @a=split /\t/,$database_site{$key2};
				$cut[8]=$a[0];  $cut[13]=$a[1];  $cut[14]=$a[2];  $cut[9]=$a[3];
			}
			my $linen=join("\t",@cut);

			next if(!exists $genelist{$cut[6]} );
			######################添加位点信息####################
			#$cut[55] = join"\|", ((&get_jiyinmiaoshu($cut[6]))[1],$cut[55]);
			my @disease_system=map{$disease{$_} || "-"}@all_gene; ### panqi 180524 
			my $diseaseasystem=join"\|",@disease_system;     #添加所有的疾病系统  
			#my @site_dase_info=&add_sitebase_info($linen);       #添加解读库信息
			#my @Disease_Namein=(&get_jiyinmiaoshu($cut[6]))[7,1];    #更新筛选的疾病系统

			my $cont_ditacion_term=$control_ditaction_ratio{join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-',1256); #NOVO1256样本检出率
			my $case_ditacion_term=$case_ditacion_ratio{join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-','-'); 
			#高血压样本检出率（实时更新）
			#say "订单${sample_id}中${key}位点未在高血压位点检出率文件（路径为：${control_detection_path}）中存在，请核查${sample_id}样本是否执行过mergevcf.sh程序。如未执行，请执行此程序；如有执行，请核查高血压config文件（路径：${hypers_cfg_addsex_path}）最新的文件中是否存在此订单，而高血压检出率文件中未找到此订单。如有此种情况，请核查mergevcf.sh的正确，并在高血压config文件去掉此订单，并重新运行mergevcf.sh程序" if ($case_ditacion_term eq "-\t-\t-\t-");
        
			my @all_sample_case=split/,/,((split/\t/,$case_ditacion_term)[2]);
			my @all_sample  =map{$doctor_patient_info{$_}->[0] // '-'}@all_sample_case;
			my @all_dotcor  =map{$doctor_patient_info{$_}->[1] // '-'}@all_sample_case;
			my @all_hospital=map{$doctor_patient_info{$_}->[2] // '-'}@all_sample_case;

			my @all_sample_gender=map{$sampleid_sex{$_} // '-'}@all_sample_case;
			my $patient_doctor_hos=(@{$doctor_patient_info{$sample_id}}) ? $doctor_patient_info{$sample_id} : [("-")x3]; #医生患者信息
			my $patient_sex_inform=$sampleid_sex{$sample_id} // "-";  #性别信息
			my $case_control_rates=
			   ((split/\t/,$case_ditacion_term)[0]=~/[0-9]+\.[0-9]*/i and (split/\t/,$case_ditacion_term)[0]>0 and
				(split/\t/,$cont_ditacion_term)[0]=~/[0-9]+\.[0-9]*/i and (split/\t/,$cont_ditacion_term)[0]>0 )
			   ?(split/\t/,$case_ditacion_term)[0]/(split/\t/,$cont_ditacion_term)[0] : '-' ;              #高血压样本检出率 和 NOVO1256样本检出率 比例

			my $hyper_last_af_case=$hyper_last_af{join"\t",@cut[0,1,3,4]} // join"\t",(("-")x4); # NOVO1256样本AF
			my $novo_1256_af_contr=$novo_1256_af{join"\t",@cut[0,1,3,4]} // join"\t",(("-")x4);  # 高血压样本AF（未实时更新，需要增加）
			#say "订单${sample_id}中${key}位点未在高血压阳性样本AF文件（$HYper_new_AF），文件未实时更新。目前最新更新市建委20170416" if ($novo_1256_af_contr eq "-\t-\t-\t-");
			
			my @all_site_dase_info;
			@all_site_dase_info=&right_site($linen,$sample_id);   ### 一个基因对多种疾病 by panqi 190311
			for(my $k=0;$k<@all_site_dase_info;$k++)                 
			{
				my $site_dase_info=$all_site_dase_info[$k];                    #添加解读库信息
				#my $site_dase_info_14=(split "\t", $all_site_dase_info[$k])[14];
				my $site_dase_info_6=(split "\t", $all_site_dase_info[$k])[6];
				my @final_out=(
					$sample_id,$linen,$diseaseasystem,$site_dase_info,
					@$patient_doctor_hos,$patient_sex_inform,
					(join",",@all_sample),(join",",@all_dotcor),(join",",@all_hospital),(join",",@all_sample_gender),
					$cont_ditacion_term,$case_ditacion_term,$case_control_rates,$novo_1256_af_contr,$hyper_last_af_case
				);
			
				my $all_output=join"\t",@final_out; #print "##$sample_id,$_\t####$final_out[4]\t$final_out[7]\t$final_out[8]\t$final_out[11]\t$final_out[12]\t$final_out[13]\n";# if(!$final_out[4]);
				$all_output=~s/多发性内分泌腺瘤综合征2型/多发性内分泌腺瘤病2型/g;
				$all_output=~s/多发性内分泌腺瘤2型/多发性内分泌腺瘤病2型/g;
				
				say ALL $all_output;
        
				##############过滤准备###############################
				#my @dp4all=(split/=|,/,$cut[26])[1..4];
				#my $dp4sum=sum @dp4all;
				#my $dp4_ratio=$cut[27];
				my $mu_type="$cut[8]";
			
				my $ratio1000=0;

				#say "$cut[19]\t$cut[20]\t$cut[21]\t$cut[22]\n";
				for(@cut[19,21,22]){($_ eq "." || $_ <= 0.01) && ++$ratio1000};
				#if($cut[16] eq "." || $cut[16] < 0.01)	{++$ratio1000;}
				#if($cut[17] eq "." || $cut[17] < 0.01)  {++$ratio1000;}

				my $sample_detected_num_1256= (split/\t/,$cont_ditacion_term)[1];
				my $sample_detected_num_case= (split/\t/,$case_ditacion_term)[1];
				my @protoin_predict=split/;/,$cut[34];
				my $protoin_predict_num=0;
				my ($SIFT_Polyphen2,$CADD)=(@protoin_predict>=4) ? ( (join"",@protoin_predict[0..3]),$protoin_predict[-1] ) : ($cut[34],$cut[34]); ##changed by sunxq 20171122
				
				##############过滤###############################
				my @filter_type=();
				my @filter_type2=();
				my $hgmd_flag=0;
				my $Disease_Namein_flag=0;
				if($cut[45]<4){push @filter_type, 'DP<4';push @filter_type2, 'DP<4';}
				#($dp4_ratio<0.1) && (push @filter_type, 'DP4_ratio<0.1');
				if($ratio1000<3){push @filter_type, 'AF1000_4>0.01';push @filter_type2, 'AF1000_3>0.01';}
				if($site_dase_info_6 eq '-'){push @filter_type, "No_disease";$Disease_Namein_flag=1;}
				if($patient_sex_inform=~/男/ and $site_dase_info_6=~/妊娠高血压综合/)
				{
					push @filter_type, 'Man_pregnancy_hypertension';
					push @filter_type2,'Man_pregnancy_hypertension';
				}
				#($sample_detected_num_case ne '-' and $sample_detected_num_case>100) && (push @filter_type, "Number_case>100");        

				#($mu_type eq "synonymous SNV") && (push @filter_type, 'synonymousSNV');	
				#($sample_detected_num_1256 ne '-' and $sample_detected_num_1256>3)   && (push @filter_type, "Number_control>3");
				#(@protoin_predict>=5 and $SIFT_Polyphen2!~/deleterious|damaging|\-/i and $CADD ne '-' and $CADD ne '.' and $CADD<10) && (push @filter_type, "protoin_predict_filter");
	
				$protoin_predict_num=0;
				for(@protoin_predict)
				{
					if($_=~/tolerated|benign/i)
					{
						$protoin_predict_num++;
					}
					if($_=~/deleterious|damaging/i){$protoin_predict_num=-1;} #wph add 190115	
				}
				if($protoin_predict_num>=3){push @filter_type, "protoin_predict_filter";}
			
				my $flag_dm;
				my $flag_dm_type;
				my ($FLAG_DM_62,$FLAG_DM_64,$FLAG_DM_T62,$FLAG_DM_T64)=(0,0,0,0);
				my @flag_dm_array;
				for(my $i=62;$i<66;$i++)
				{
					next if($i==63 || $i==65);
					next if($cut[$i] eq "-");
					@flag_dm_array=split /\|/,$cut[$i];
					for(@flag_dm_array)
					{
						$flag_dm_type=(split /:/,$_)[0];
						next if($flag_dm_type ne "true_type");
						$FLAG_DM_T62=1 if($i==62);
						$FLAG_DM_T64=1 if($i==64);
				
						$flag_dm=(split /:/,$_)[3];	
						if($flag_dm!~/DM|pathogenic|significance/i && $flag_dm ne "not provided" && $flag_dm ne "other" )
						{
							$FLAG_DM_62=1 if($i==62);
							$FLAG_DM_64=1 if($i==64);
						}
						#if($flag_dm=~/DM|pathogenic/i && $flag_dm!~/conflict/i )
						#{
						#	$FLAG_DM_62=2 if($i==62);
						#	$FLAG_DM_64=2 if($i==64);
						#}
					}
				}
			
				if($mu_type eq "nonsynonymous SNV" || $mu_type=~/splicing/ )
				{
					if( ($FLAG_DM_T62==1 && $FLAG_DM_T64==1 && $FLAG_DM_62==1 && $FLAG_DM_64==1) || ($FLAG_DM_T62==1 && $FLAG_DM_T64==0 && $FLAG_DM_62==1 ) || ($FLAG_DM_T62==0 && $FLAG_DM_T64==1 && $FLAG_DM_64==1) )
					{
						$hgmd_flag=1;
						push @filter_type2, "non_pathogenic_DB";
					}elsif($FLAG_DM_T62==1 || $FLAG_DM_T64==1){
						$hgmd_flag=1;
					}
				}

				if( $mu_type eq "synonymous SNV" || $cut[5] eq 'intronic' || $cut[5]=~/UTR/ )
				{
					if( ($FLAG_DM_T62==1 && $FLAG_DM_62!=1) || ($FLAG_DM_T64==1 && $FLAG_DM_64!=1) )
					{
						$hgmd_flag=1;
					}else{
						$hgmd_flag=1;
						push @filter_type2, "synonymous_SNV" if($mu_type eq "synonymous SNV");
						push @filter_type2, "intronic_syn" if($cut[5] eq 'intronic');
						push @filter_type2, "UTR" if($cut[5]=~/UTR/);
					}
				}
				#if(($cut[5] eq 'intronic' || $cut[8] eq 'synonymous SNV' ) && $FLAG_DM_T62==0 && $FLAG_DM_T64==0 ){ #sxq changed 190508
                                #	push @filter_type, "intronic_syn";
                        	#}	
				#$hgmd_flag=1 if( $cut[43]=~/true/i && $cut[43]=~/DM/i );  ### panqi 190124
			
				my $flag=(@filter_type) ? (join",",@filter_type) : '?';
				my $filter_type_size=@filter_type;
				my $filter_type2_size=@filter_type2;
				#if (($cut[52] ne "-" and $flag!~/Man_pregnancy_hypertension/i ) or ($hgmd_flag==1 && $ratio1000>2 && $Disease_Namein_flag==0 && @filter_type<1)  or @filter_type<1){  ### panqi 181017 add ratio1000 判断
				if ( ($cut[52] ne "-" and $flag!~/Man_pregnancy_hypertension/i ) or ($hgmd_flag==1 && $filter_type2_size==0) or ($hgmd_flag==0 && $filter_type_size==0) ){
					++$sample_number{$sample_id};
					if(exists $hash_base{$key2})
					{
						$flag="done";
					#	$flag.="_NM" if($hash_base2{$key2} eq "-");
					}
					else{	$flag="todo"; }
					#$flag=(exists $hash_base{$key2}) ? "done" :"todo";
				
					if(exists $hash_base_fubiao{$key2} && !exists $hash_base{$key2})
					{
						$flag="done_fubiao";
					#	$flag.="_NM" if($hash_base_fubiao2{$key2} eq "-");
					}
					#$flag="done_fubiao" if(exists $hash_base_fubiao{$key2} && !exists $hash_base{$key2});
				
					$flag="todo_oldDB" if(exists $right_old{join"\t",@cut[0,1,3,4]}  && !exists $hash_base{$key2} && !exists $hash_base_fubiao{$key2});  #######sunxq 20171205
					if($hgmd_flag==1){
						if(($cut[52] ne "-" and $flag!~/Man_pregnancy_hypertension/i) or @filter_type<1){
						}else{
					#		$flag.="_DM";
						}
					}
					say RIGHT "$flag\t$all_output";
					#if ($flag eq "base_to_read"){
					#say READ  join"\t",("-\t-\t-",(split/\t/,$case_ditacion_term)[2],$_,$diseaseasystem,@site_dase_info) ;
					#}
				}else{
					if($hgmd_flag==1){ $flag=(@filter_type2) ? (join",",@filter_type2) : '?';}
					if($hgmd_flag==0){ $flag=(@filter_type) ? (join",",@filter_type) : '?';}
					if(exists $hash_base{$key2}){$flag.=";exist_in_readbase";}
					if(exists $hash_base_fubiao{$key2} && !exists $hash_base{$key2}){$flag.=";exist_in_FubiaoReadbase";}
					if($Disease_Namein_flag==0 || $cut[43]=~/true_positive/i){
						say WRONG "$flag\t$all_output";                    #### panqi 180406 ???
					}
				}
			} 
		}
		close IN;
	}
	foreach my $idsamle (keys%sample_number){
		($sample_number{$idsamle}<=2) && (say RIGHT join"\t",('Negative_sample',"$idsamle",('-')x103));  #### panqi 180406 ???
	}
	say "step 1执行完成！！";
}
close FI;
close ALL;
close RIGHT;
close WRONG;


sub sub_step_2{
#&sub_hash_base;
#&sub_hash_deletion;
my %sampledata=&read_samdatabase;
my %coresitedb;
my %sample_flag;

open RIGHT,   ">>$out/GXY_sample_rigth_info.txt"   or die "Cannot find the File :$out/GXY_sample_rigth_info.txt";
open WRONG,   "$out/GXY_sample_wrong_info.txt"   or die "Cannot find the File :$out/GXY_sample_wrong_info.txt";
while(<WRONG>){
	next if /^Flag/i;
	chomp;
	my @cut=split/\t/;
	my $n=@cut;
	next if( $n < 20);
        my $key;
        if($cut[45]=~/true/i){
                  $key=join("_",@cut[1..3,5,6,83]); #add end HGMD id 2016.08.23
        }else{
                if($cut[46]=~/snp/i || $cut[46] eq '-'){     # add by sxq 2017.8.15
                          $key=join("_",@cut[1..3,5,6],$cut[83]);
                }else{
                          $key=join("_",@cut[1..3,5,6],$cut[83]);
                }
        }
	print RIGHT "$_\n" if(exists $sampledata{$key}{key});     #### panqi 180406 ???
}
close RIGHT;
close WRONG;


open RIGHT,   "$out/GXY_sample_rigth_info.txt"   or die "Cannot find the File :$out/GXY_sample_rigth_info.txt";
open RIGHTFIN,">$out/GXY_sample_rigth_info_final.txt" or die $!;
#open CORESITE,">$out/GXY_sample_rigth_core_site.txt" or die $!;
open FINALPDF,">$out/GXY_sample_rigth_to_PDF_final.txt" or die $!;
say RIGHTFIN join"\t",("Flag",@titel);
#say CORESITE join"\t",@titel[1..4,6..$#titel];
while(<RIGHT>){
	next if /^##/i;
	chomp;
	my @cut=split/\t/;
	(say FINALPDF join"\t",@cut[1..70,76,71..75,55,77..89,90,-2,-1] and next) if /^Flag/i;
	my $key;
	if($cut[5] eq '.'){$cut[5]='-';}
	if($cut[45]=~/true/i){
		$key=join("_",@cut[1..3,5,6,83]); #add end HGMD id 2016.08.23
	}else{
		if($cut[46]=~/snp/i || $cut[46] eq '-'){     # add by sxq 2017.8.15
			$key=join("_",@cut[1..3,5,6],$cut[83]);
		}else{
			$key=join("_",@cut[1..3,5,6],$cut[83]);
		}
	}
	print "$key\n";
	my $linenew='';
	if(exists $sampledata{$key}{key}){ 
		$sample_flag{$cut[1]}{k}=1;  
		if($sampledata{$key}{key} =~/report/i || $sampledata{$key}{key} =~/ok/i || $sampledata{$key}{key} =~/Negative/i){
			$linenew="$sampledata{$key}{info}";
			$sample_flag{$cut[1]}{f}=1;  
		}elsif($sampledata{$key}{key} =~/del/i){
			$linenew="####$sampledata{$key}{info}";
		}elsif($sampledata{$key}{key} =~/beizhu/i){
			$linenew="###$sampledata{$key}{info}";
		}elsif($sampledata{$key}{key} =~/fubiao/i){
			$linenew="##$sampledata{$key}{info}";
		}else{
			die "请核实$key\备注为$sampledata{$key}{key}\n";
		}
        }else{
                shift @cut;
                my $slnew=join("\t",@cut);
                $linenew="N$slnew\n";
        }
        my @cutnew=split /\t/,$linenew;
        unshift @cutnew,'Moo';
	#(say FINALPDF join"\t",@cut[1..56,92,58..70,76,71..75,55,77..89] and next) if /^Flag/i;
        my $key3=join"\t",@cutnew[2,3,5,6];
	#$coresitedb{$key3}=[@cutnew[2..89,96..$#cutnew]] unless /^Negative_sample/i;
	($cutnew[78]=~/良性$/i) && (say "订单号$cutnew[2]的${key3}位点检查为$cutnew[78]位点，其主要是是由于数据挖掘组和产品经理制定的蛋白预评定阈值不一致，请及时和产品经理沟通！！");
	say RIGHTFIN join"\t",@cutnew;
	say FINALPDF join"\t",@cutnew[1..70,76,71..75,55,77..89,90,118];#wph add 解读人20181207
}
#map{say CORESITE join"\t",@$_}(values %coresitedb);
say "step 2执行完成！！";

foreach my $key(keys %sample_flag){
	if(!exists $sample_flag{$key}{f} && exists $sample_flag{$key}{k}){
		my $null="-\t" x 89;
		print FINALPDF "$key\t$null\n";
	}
}

close RIGHT;
close RIGHTFIN;
#close CORESITE;
close FINALPDF;
}

##########################子函数区域I###########################################
###########################高血压样本中，阳性位点的样本检出率####################
sub sub_case_ditacion_ratio{
    my $detection_ratio_prefix="snp_indel_detection_out";
    my $latest_detection_file=judge_lastest_sitedb($case_detection,$detection_ratio_prefix);#查找最新的样本检出率文件
say '*'x25;
say "正在读取最新高血压检出率文件： 
$latest_detection_file";
    open CASEDI,"$latest_detection_file" || die "can't open the $latest_detection_file!";
        %case_ditacion_ratio=map{chomp;(join"\t",(split/\t/,$_)[0..3]) => (join"\t",(split/\t/,$_)[4..7])}<CASEDI>;
    close CASEDI;
say '读取完成！！';
}
###################################### 阴性对照检出率 #################################################
sub sub_control_ditaction_ratio{
    my $control_detection_prefix="detection_ratio_out";
    my $latest_control_detection_ratio=judge_lastest_sitedb($control_detection_path,$control_detection_prefix);
say '*'x25;
say "正在读取最新高血压1256NOVO1000阴性样本检出率文件： 
$latest_control_detection_ratio";
    open CONTR,"$latest_control_detection_ratio" || die "can't open the $latest_control_detection_ratio!";
        %control_ditaction_ratio=map{chomp;(join"\t",(split/\t/,$_)[0..3]) => (join"\t",(split/\t/,$_)[4..7])}<CONTR>;
    close CONTR;
say '读取完成！！';
}
############################# 高血压样本信息路径 ################################################## 
sub sub_doctor_patient_info{
    my $hyper_sample_info_prefix="hypers_sample_info";
    my $latest_hyper_sample_info=judge_lastest_sitedb($hyper_sample_info_path,$hyper_sample_info_prefix);
say '*'x25;
say "正在读取最新高血压医生患者信息单文件： 
$latest_hyper_sample_info";
    open DOCT,"$latest_hyper_sample_info" || die "can't open the $latest_hyper_sample_info!" ;
        %doctor_patient_info=map{chomp;(split/\t/,$_)[0] => [(split/\t/,$_)[1..3]]}<DOCT>;
    close DOCT;
say '读取完成！！';
}
########################## 疾病的指导意见 #########################################################
sub sub_gene_disease{
    my $gene_desease_prefix="gene_desease";
    my $latest_gene_desease=judge_lastest_sitedb($gene_desease_path,$gene_desease_prefix);
    
say '*'x25;
say "正在读取最新基因——疾病文件： 
$latest_gene_desease";
    open GEDI,$latest_gene_desease or die "can't open the file:$latest_gene_desease";
        map{chomp; $gene_disease{(split/\t/,$_)[0]}{(split/\t/,$_)[1]}=join"\t",(split/\t/,$_)[1..7];}<GEDI>;
	#map{chomp; $genelist{(split/\t/,$_)[0]}=1; }$line;
    close GEDI;
say '读取完成！！';
}
###########################  gene list  ###################################
sub sub_gene{
	my $gene_desease_prefix="gene_desease";
	my $latest_gene_desease=judge_lastest_sitedb($gene_desease_path,$gene_desease_prefix);

say '*'x25;
say "正在读取最新基因——疾病文件：
$latest_gene_desease";
	open GEDI,$latest_gene_desease or die "can't open the file:$latest_gene_desease";
		map{chomp; $genelist{(split/\t/,$_)[0]}=1; }<GEDI>;
	close GEDI;
say '读取完成！！';	
}
########################## 基因疾病系统 #######################################################
sub sub_disease{
    open DISEA,"$disease_file" || die "can't open the file:$disease_file";
        %disease=map{chomp;(split/\t/,$_)[0]=>(split/\t/,$_)[1];}<DISEA>;
    close DISEA;
}

sub sub_disease_ratio{
	open DISEA,"$disease_ratio_file" || die "can't open the file:$disease_ratio_file";
		%disease_ratio=map{chomp;(split/\t/,$_)[0]=>(split/\t/,$_)[1];}<DISEA>;
	close DISEA;
}
########################## 高血压位点库文件 ###################################################

sub sub_hash_base{
    my $latest_gxybase_prefix=".right.txt";
    my $latest_gxybase=&get_latest_data($latest_gxybase_path,$latest_gxybase_prefix);
	
    my $key;
    my $info;
say '*'x25;
say "正在读取最新高血压位点库文件： 
$latest_gxybase";
    open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
    while(<DATA>){
	$key=join"\t",(split/\t/,$_)[3,4,6,7];
	$info=join"\t",(split/\t/,$_)[78..92];
	push @{ $hash_base{$key} }, $info;
	
	$info=join"\t",(split/\t/,$_)[11..14];
	$database_site{$key}=$info;
	$hash_base2{$key}=(split/\t/,$_)[91];
	
	$key=join"\t",(split/\t/,$_)[2,3,4,6,7];
	$hash_base_2{$key}=1;
    }
    close DATA;

say '读取完成！！';
}

##################高血压附表位点库########3
sub sub_hash_fubbase{
	my $latest_gxybase_prefix=".right.txt";
	my $latest_gxybase=&get_latest_data($latest_gxybase_path_fubiao,$latest_gxybase_prefix);

	my $ID;
say '*'x25;
say "正在读取最新高血压附表位点库文件：
$latest_gxybase";
	open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
	while(<DATA>){
		chomp;
		my $line=$_;
		next unless($line);
		#next if($line=~/^$/);
		$ID=(split/\t/,$line)[2];
		my $key=join"\t",(split/\t/,$line)[3,4,6,7];
		my $info=join"\t",(split/\t/,$line)[78..90];
		push @{ $hash_base_fubiao{$key} }, $info;
		#$hash_base_fubiao2{$key}=(split/\t/,$line)[91];
	}
	close DATA;	
}
############################高血压出报告位点库文件#######

sub read_samdatabase{
	my $latest_gxybase_prefix=".right.txt";
	my $latest_gxybase=&get_latest_data($latest_gxyreportbase_path,$latest_gxybase_prefix);
say '*'x25;
say "正在读取最新高血压出报告位点库文件： 
$latest_gxybase";
	open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase"; 
	my %hash; 
	while(<DATA>){
		chomp;
		my @tmp=split/\t/,$_;
		my $n=@tmp;
		$n=$n-1;
		if($tmp[6] eq '.'){$tmp[6]='-';}
		if($tmp[7] eq '.'){$tmp[7]='-';}
		my $keysam=join("_",@tmp[2..4,6,7,84]);
		my $info=join("\t",@tmp[2..$n]);
		$hash{$keysam}{info}=$info;
		$hash{$keysam}{key}=$tmp[0]; 
	}
	close  DATA;
	return %hash; 
	say '读取完成！！';
}

##########################################################
sub sub_hash_oldbase{
	my $key;
	my $info;

say '*'x25;
say "正在读旧版取高血压位点库文件：
$latest_gxybase_old_path";
	open DATA,$latest_gxybase_old_path or die "can't open the file:$latest_gxybase_old_path";
	while(<DATA>)
	{
		$key=join"\t",(split/\t/,$_)[4,5,7,8];
		$info=join"\t",(split/\t/,$_)[67..79];
		$right_old{$key}=$info;
	}
	close DATA;
}

###########################################################
sub sub_hash_deletion{
    my $latest_delbase_prefix=".right.txt";
    my $latest_delbase=&get_latest_data($latest_deletion_path,$latest_delbase_prefix);
say '*'x25;
say "正在读取最新待删除位点库:
$latest_delbase";
    open DATA,$latest_delbase or die "can't open the file:$latest_delbase";
	%hash_del=map{chomp; (join"\t",(split/\t/,$_)[1,2,3,5,6]) => 1  ;}<DATA>;		         
    close DATA;
say '读取完成！！';
}
#############################高血压样本性别config文件#######################
sub sub_sampleid_sex{
    my $hypers_cfg_addsex_prefix="hypers_cfg_addsex.xls";
    my $latest_hypers_cfg_addsex=judge_lastest_sitedb($hypers_cfg_addsex_path,$hypers_cfg_addsex_prefix);
say '*'x25;
say "正在读取最新高血压样本性别config文件： 
$latest_hypers_cfg_addsex";
    open CONF,"$latest_hypers_cfg_addsex" || die "can't open the $latest_hypers_cfg_addsex!";
        %sampleid_sex=map{chomp; (split/\t/,$_)[0] => (split/\t/,$_)[5];}<CONF>;
    close CONF;
say '读取完成！！';
}
#############################高血压NOVO1256AF#######################
sub sub_novo_1256_af{
say '*'x25;
say "正在读取最新高血压样NOVO_1256_AF文件： 
$NOVO_1256_AF";
    open AFNOVO,"$NOVO_1256_AF" || die "can't open the $NOVO_1256_AF!";
        %novo_1256_af=map{chomp; (join"\t",(split/\t/,$_)[0,1,3,4])=>(join"\t",(split/\t/,$_)[5..8]);}<AFNOVO>;
    close AFNOVO;
say '读取完成！！';
}
#############################高血压case_last_AF#######################
sub  sub_hyper_last_af{
say '*'x25;
say "正在读取最新高血压样HYper_new_AF文件： 
$HYper_new_AF";
    open AFHYPE,"$HYper_new_AF" || die "can't open the $HYper_new_AF!";
    %hyper_last_af=map{chomp; (join"\t",(split/\t/,$_)[0,1,3,4])=>(join"\t",(split/\t/,$_)[5..8]);}<AFHYPE>;
    close AFHYPE;
say '读取完成！！';
}

##########################子函数区域II############################
sub get_latest_data{
    my ($path,$base_name)=@_;
    my @all_file_name=grep{chomp; -d "$path\/$_" }(`ls $path`);
    my $max_file_name=max (@all_file_name);  
    #my $latest_file=glob("$path/$max_file_name/*$base_name*"); 
    my $latest_file=`ls $path/$max_file_name/*$base_name`;
    return $latest_file;
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

=head
sub get_jiyinmiaoshu{
    my $gene=shift;
    my $cha_gene=$gene=~s/\(.*\)//gr; 
    my @all_gene=split/;/,$cha_gene;

    my @chose_gene;
    foreach my $eachgene (@all_gene){
        if (exists $gene_disease{$eachgene}){
            foreach my $eachdisease (keys %{$gene_disease{$eachgene}}){
                my @each_term=($eachgene,(split/\t/,$gene_disease{$eachgene}{$eachdisease}));
                push @chose_gene,\@each_term;    
            }
        }
    }
    
    @chose_gene=([$cha_gene,("-")x7]) if (!@chose_gene);

    my  @chose_gene_new;
    for (my $i=0;$i<@chose_gene;$i++){
        for (my $j=0;$j<@{$chose_gene[$i]};$j++){
            $chose_gene_new[$j][$i]=$chose_gene[$i][$j];
        }
    }
    my @all_term=map{join"\|",@$_}@chose_gene_new;
    return @all_term;
}
=cut
 
sub get_jiyinmiaoshu{
    my $gene=shift;
    my $cha_gene=$gene=~s/\(.*\)//gr; 
    my @all_gene=split/;/,$cha_gene;
	my %hash=();
	my @all_term=();

	my @chose_gene;
	foreach my $eachgene (@all_gene){
		if (exists $gene_disease{$eachgene}){
			foreach my $eachdisease (keys %{$gene_disease{$eachgene}}){
				push @{ $hash{$eachgene} }, $eachgene."\t".$gene_disease{$eachgene}{$eachdisease}; 
				push @all_term, $eachgene."\t".$gene_disease{$eachgene}{$eachdisease};
			}

		}
		else{
			push @{ $hash{$eachgene} }, "$eachgene\t-\t-\t-\t-\t-\t-\t-";
			push @all_term, "$eachgene\t-\t-\t-\t-\t-\t-\t-";
		}
	}
	return @all_term;   
}
 
#	foreach my $eachgene (@all_gene)
#	{
#		my @arr=@{$hash{$eachgene}};
#		for( my $i=0; $i<@arr; $i++ )
#		{
#			push @all_term, $arr[$i];
#		}
#	}


sub right_site{
	my $site_info=shift;
	my $sample = shift;
	next if /^CHROM/;
	my @cut=split/\t/, $site_info;
	my @tmp;
	my @add_info_all=();
	my %flag_add_info=();
	my ($gene,$jibing);
            
	my @add_info=(("-")x16);
	my $key;
	my $line="";
	if($cut[43]=~/true/i){
		$key=join("\t",@cut[0,1,3,4]); #add end HGMD id 2016.08.23
	}else{
		if($cut[44]=~/snp/i){     # add by sxq 2017.8.15
			$key=join("\t",@cut[0,1,3,4]);
		}else{
			$key=join("\t",@cut[0,1,3,4]);
		}
	}
    
	if( exists $hash_base{$key} )
	{
		@tmp=@{ $hash_base{$key} };
		for(my $k=0;$k<@tmp;$k++)
		{
			$gene=(split/\t/,$tmp[$k])[8];
			$jibing=(split/\t/,$tmp[$k])[6];
			@add_info[0..14] = (split/\t/,$tmp[$k])[0..14];
			$flag_add_info{$gene}{$jibing} = 1;
			###@add_info[8,6,2,7,11,3,12]= (&get_jiyinmiaoshu($cut[6]))[0,1,2,3,4,5,6];
			$add_info[15]=$add_info[14];
			$add_info[14]=$add_info[13];
			$add_info[13]=$cut[12]; #CDS    PEP
			if($cut[30]=~/het/ and $add_info[7]=~/AR/ and $add_info[7]!~/AD/){
				$add_info[4]="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，".$add_info[4];
			}#wph add 190104
			$line=join("\t",@add_info[0..15]);
			push @add_info_all, $line;
		}
	}
	if( !exists $hash_base{$key} && exists $hash_base_fubiao{$key}){
		@tmp=@{ $hash_base_fubiao{$key} };
		for(my $k=0;$k<@tmp;$k++)
		{
			$gene=(split/\t/,$tmp[$k])[8];
			$jibing=(split/\t/,$tmp[$k])[6];
			@add_info[0..12] = (split/\t/,$tmp[$k])[0..12];
			$flag_add_info{$gene}{$jibing} = 1;

			my @tmp2=&get_jiyinmiaoshu($gene);
			for(my $k2=0;$k2<@tmp2;$k2++)
			{
				$add_info[14]=(split/\t/,$tmp2[$k2])[7] if($jibing eq (split/\t/,$tmp2[$k2])[1]);
				$add_info[15]=(split/\t/,$tmp2[$k2])[1] if($jibing eq (split/\t/,$tmp2[$k2])[1]);
			}
			$add_info[13]=$cut[12]; #CDS    PEP
			if($cut[30]=~/het/ and $add_info[7]=~/AR/ and $add_info[7]!~/AD/){
				$add_info[4]="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，".$add_info[4];
			}#wph add 190104
			$line=join("\t",@add_info[0..15]);
			push @add_info_all, $line;
		}	
	}
	
	if (exists $gene_disease{$cut[6]})
	{
		foreach my $eachdisease (keys %{$gene_disease{$cut[6]}})
		{
			if( !exists $flag_add_info{$cut[6]}{$eachdisease})
			{
				@add_info=&add_sitebase_info($site_info, $eachdisease, $sample);
				$line=join("\t",@add_info[0..15]);
				push @add_info_all, $line;
			}
		}
	}
        
	
	return @add_info_all;
}
close IN;

sub add_sitebase_info{
    my $site_info=shift;
    my $disease=shift;
    my $sample=shift;
    next if /^CHROM/;
    my @cut=split /\t/, $site_info;
	my @tmp;
	my @add_info_all=();
            
    my @add_info=(("-")x16);
    my $key;
    if($cut[43]=~/true/i){
		$key=join("\t",@cut[0,1,3,4]); #add end HGMD id 2016.08.23
    }else{
		if($cut[44]=~/snp/i){     # add by sxq 2017.8.15
			$key=join("\t",@cut[0,1,3,4]);
		}else{
			$key=join("\t",@cut[0,1,3,4]);
		}
    }

	if($cut[0]==5 && $cut[1]==136969729)
	{
		$cut[1]=136969729;
	}
    
        my $mutation_type=
            ($cut[8]=~/stopgain/i)      ? "终止"  :
            ($cut[8]=~/splicing/i)      ? "剪切":	
            ($cut[8]=~/nonframeshift/i) ? "非移码":
            ($cut[8]=~/frameshift/i)    ? "移码"  :
            ($cut[8]=~/nonsynonymous/i) ? "错义"  :
            ($cut[8]=~/synonymous/i)    ? "同义"  : "?";
        my $detection_or_not=
            ($cut[19] eq "." and $cut[20] eq ".") ? "均未检出" :'均低于1%';
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
               $ratio_detection.="千人项目东亚人群检出频率为$ratio_chinese，ExAC数据库东亚人群检出频率为$ratio_ExAC，gnomAD数据库东亚人群检出频率为$ratio_gnomAD，本地数据库检出频率为$ratio_novo。";
        }else{
              if($ratio_chinese ne '未收录'){$ratio_detection.="千人项目东亚人群检出频率为$ratio_chinese，";}
              if($ratio_ExAC ne '未收录'){$ratio_detection.="ExAC数据库东亚人群检出频率为$ratio_ExAC，";}
              if($ratio_gnomAD ne '未收录'){$ratio_detection.="gnomAD数据库东亚人群检出频率为$ratio_gnomAD，";}
              if($ratio_novo ne '未收录'){$ratio_detection.="本地数据库检出频率为$ratio_novo，";}
              my $l=join("、",@detection_ratio_n);
              $ratio_detection.="$l\未收录。";
        }

        my $software_pred='';
        my $splicing_pred='';
        #    ($cut[34] eq '.' or $cut[34] =~ /-;-;-;-;/) ? "" : "，软件预测该突变对基因或基因产物有害";
        my @softwares=split /;/,$cut[34];
        my @ssss;
        my $sss_n=0;
        if($softwares[0] =~/Deleter/i || $softwares[0] =~/dama/i){push @ssss,'SIFT';$sss_n++;}  #SIFT、Polyphen2_HVAR、Polyphen2_HDIV、M-CAP
        if($softwares[1] =~/Deleter/i || $softwares[1] =~/dama/i){push @ssss,'Polyphen2_HVAR';$sss_n++;}
        if($softwares[2] =~/Deleter/i || $softwares[2] =~/dama/i){push @ssss,'Polyphen2_HDIV';$sss_n++;}
        if($softwares[3] =~/Deleter/i || $softwares[3] =~/dama/i){push @ssss,'M-CAP';$sss_n++;}
        my $sssware=join("、",@ssss);
        if($sss_n>0){
        	$software_pred="$sssware\软件预测该变异对基因或基因产物有害。";
        }
		if($cut[8]=~/splicing/i)
		{
			if( ($cut[28] ne ".") and ($cut[28] > 0.6) and ($cut[29] ne ".") and ($cut[29] > 0.6) and ($cut[32] ne ".") and abs($cut[32]) >=2)
			{
				$splicing_pred.="dbscSNV1.1、spidex软件预测该变异对基因或基因产物有害。";
			}
			elsif( ($cut[28] ne ".") and ($cut[28] > 0.6) and ($cut[29] ne ".") and ($cut[29] > 0.6) )
			{
				$splicing_pred .="dbscSNV1.1软件预测该变异对基因或基因产物有害。";
			}
			elsif( ($cut[32] ne ".") and abs($cut[32])>=2)
			{
				$splicing_pred.="spidex软件预测该变异对基因或基因产物有害。";
			}
			else
			{
				$splicing_pred.="";
			}
		}
		my $wenxianzhichi=($cut[62] eq "-" and $cut[63] eq "-" and $cut[64] eq "-" and $cut[65] eq "-" ) ? "由于暂无相关研究，" : "由于研究较少，";
		#($cut[58]=~/\S/ and (length $cut[58]) >10 and $cut[59]=~/http/) ? "由于研究较少，" : "由于暂无相关研究，";
        
		my $muty=""; #wph changed 190104 "$cut[14]${mutation_type}";
		if($cut[5]=~/splicing/){
			$muty="位于剪切区";
		}else{
			$muty="会导致$add_info[10]${mutation_type}变异";
		}
		$muty=~s/^-//;
		if(!exists $right_old{$key}){	
			@tmp=&get_jiyinmiaoshu($cut[6]);
			for(my $k=0;$k<@tmp;$k++)
			{
				next if( (split/\t/,$tmp[$k])[1] ne $disease);
				@add_info[0,1]=($cut[43]=~/true/i) ? ("待解读","生信") : ("待解读","生信") ; #ACMG条目,#致病性结论
				@add_info[8,6,2,7,11,3,12,14,15]= (split/\t/,$tmp[$k])[0..7,1];  #致病基因 疾病名称 疾病描述 遗传模式 临床指导建议 基因描述    临床表现 疾病大类
				@add_info[9,10,13]=@cut[13,14,12]; #CDS    PEP
				$add_info[10]=&pep_trans($cut[14]); #wph add 190122
				my @acmg=&acmg_auto_judge($site_info,$disease,$add_info[7],$sample);
				$add_info[0]=join(",",@acmg[0..$#acmg]) if(@acmg > 0);
				$add_info[0]="无" if(@acmg == 0);
				if($cut[5]=~/splicing/){
					$muty="位于剪切区";
				}else{
					$muty="会导致$add_info[10]${mutation_type}变异";
				}
				my $arhet="";#wph changed 190104
				if($cut[30]=~/het/ and $add_info[7]=~/AR/ and $add_info[7]!~/AD/){
					$arhet="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，";
				}
				$add_info[4]="本筛查检测出$cut[6]基因的变异，变异位点为$cut[13]，查询ClinVar等公共数据库显示，$cut[13]变异${muty}。$ratio_detection"."${software_pred}${splicing_pred}${wenxianzhichi}行ACMG标准，判定该变异的致病性为临床意义未明。${arhet}具体情况请结合临床相关资料综合判断。";#位点描述
				#该突变在千人基因组和ExAC数据库中的基因频率${detection_or_not}"."${software_pred}。${wenxianzhichi}行ACMG标准，判定该突变的致病性为临床意义未明。具体情况请结合临床相关资料综合判断。";#位点描述
				$add_info[5]="-"; #参考文献
			}
		}

		if(exists $right_old{$key}){                                     ### panqi 20181025
			@tmp=&get_jiyinmiaoshu($cut[6]);
			my $site_info_tmp;
			
			for(my $k=0;$k<@tmp;$k++)
			{
				next if( (split/\t/,$tmp[$k])[1] ne $disease);
				@add_info=split/\t/,$right_old{$key};
				$site_info_tmp=$add_info[4];
				@add_info[8,6,2,7,11,3,12,14,15]= (split/\t/,$tmp[$k])[0..7,1];
				@add_info[9,10,13]=@cut[13,14,12]; #CDS    PEP 
				$add_info[10]=&pep_trans($cut[14]); #wph add 190122
				if($cut[5]=~/splicing/){
					$muty="位于剪切区";
				}else{
					$muty="会导致$add_info[10]${mutation_type}变异";
				}
				my $arhet="";#wph changed 190104
				if($cut[30]=~/het/ and $add_info[7]=~/AR/ and $add_info[7]!~/AD/){
					$arhet="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，";
				}
				$add_info[4] = "本筛查检测出$cut[6]基因的变异，变异位点为$cut[13]，查询ClinVar等公共数据库显示，$cut[13]变异${muty}。${ratio_detection}"."${software_pred}${splicing_pred}${wenxianzhichi}行ACMG标准，判定该变异的致病性为临床意义未明。${arhet}具体情况请结合临床相关资料综合判断。"."<-->$site_info_tmp"; #位点描述
			}
		}
	
	return @add_info;
}
close IN;

sub acmg_auto_judge{                    #wph add 190312       
	my $site_info = shift;
	my $disease   = shift;
	my $adar      = shift;
	my $sample    = shift;

        my @cut=split/\t/,$site_info;
        my @acmg=();
	my ($pep1,$pep2,$pep3)=("","","");
	my @site;

        if($cut[33] == 1 && $cut[8]=~/^nonsynonymous/ ){
                push @acmg,"PP3";
        }
        if($cut[34]!~/deleterious|damaging|-/i && $cut[8]=~/^nonsynonymous/ ){
                push @acmg,"BP4";
        }

	my ($pep3_num1,$pep3_num2,$pep3_rate,$pep3_flag)=(0,0,0.0001,0);

        my $ref_alt=$cut[3]."_".$cut[4];
        my $repeat_site = $cut[15]."_".$cut[16];
        if( $cut[8]=~/^stoploss/ || ( $repeat_site eq "\._\." && $cut[8]=~/^nonframeshift/) ){
                push @acmg,"PM4";
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
		if($flag_BP7 < 3 && $cut[8]=~/^synonymous/)
		{
			push @acmg,"BP7";
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
			next if( (split /:/,$_)[0] !~ /true_site/i );
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
	if( $cut[58] ne "-" )
	{
		@site=split /\|/,$cut[58];
		$pep3_num2=@site;
		for(@site)
		{
			$pep3=(split /:/,$_)[5];
			$pep3_num1++ if($pep3 eq "DM" || $pep3 eq "DM?" );
		}
		$pep3_rate=sprintf("%.2f", $pep3_num1/$pep3_num2);
		#$pep3_rate=$pep3_num1/$pep3_num2;
		if($pep3_num2 >= 7 and $pep3_rate>=0.6){$pep3_flag=1;}
	}

	$pep3_num1=0;
	if($cut[59] ne "-" )
	{
		@site=split /\|/,$cut[59];
		$pep3_num2=@site;
		for(@site)
		{
			$pep3=(split /:/,$_)[2];
			$pep3_num1++ if($pep3=~/Pathogenic/ || $pep3=~/Likely pathogenic/ );
		}
		$pep3_rate=sprintf("%.2f", $pep3_num1/$pep3_num2);
		#$pep3_rate=$pep3_num1/$pep3_num2;
		$pep3_flag=1 if($pep3_num2>=7 && $pep3_rate>=0.6);
	}
=cut

	$pep3_flag=0;
	if( $cut[59] ne "-" && $cut[8]=~/^nonsynonymous/ )
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

=head
	if( $cut[8]=~/^stopgain/ || $cut[8]=~/^frameshift/ )
	{
		if( $cut[56] ne "-" )
		{
			@site=split /\|/,$cut[56];
			$pep3_num2=@site;
			for(@site)
			{
				$pep3=(split /:/,$_)[5];
				$pep3_num1++ if($pep3 eq "DM" || $pep3 eq "DM?" );
			}
			$pep3_rate=sprintf("%.2f", $pep3_num1/$pep3_num2);
			$pep3_flag=1 if($pep3_num1>0);
		}
		$pep3_num1=0;
		if( $cut[57] ne "-" )
		{
			@site=split /\|/,$cut[57];
			$pep3_num2=@site;
			for(@site)
			{
				$pep3=(split /:/,$_)[2];
				$pep3_num1++ if( $pep3=~/Pathogenic/ || $pep3=~/Likely pathogenic/ );
			}
			$pep3_rate=sprintf("%.2f", $pep3_num1/$pep3_num2);
			$pep3_flag=1 if($pep3_num1>0);
		}
		if($pep3_flag==1){push @acmg,"PVS";}	
	
	}
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

	if( ($cut[19] eq "." || $cut[19] eq "0") && ($cut[21] eq "." || $cut[21] eq "0") && $db_dep_flag20==1 && $adar=~/AD/ )
	{
		$flag_PM2=1;
	}
	
	my $cont_ditacion_term=$control_ditaction_ratio{ join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-',1256); #NOVO1256样本检出率
	my $case_ditacion_term=$case_ditacion_ratio{     join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-','-');
	if( ( ( (split/\t/,$case_ditacion_term)[1]==1 && (split/\t/,$case_ditacion_term)[2] eq $sample ) || (split/\t/,$case_ditacion_term)[1]==0 )  && (split/\t/,$cont_ditacion_term)[1]==0 )
	{
		$flag_PM2=1 if( $adar=~/AD/ && $cut[44] eq "snp");
	}
	#if($flag_PM2==1){ push @acmg,"PM2"; }
	
	my $ratio1=$cut[19];
	$ratio1=0 if($cut[19] eq "." && $db_dep_flag1 eq "Y" && $db_dep_12>=20 );
	my $ratio2=$cut[21];
	$ratio2=0 if($cut[21] eq "." && $db_dep_flag2 eq "Y" && $db_dep_22>=20 );
	my $ratio3=$cut[20];
	$ratio3=0 if($cut[20] eq ".");
	my $ratio_max = 1;
	if($ratio1 ne "." && $ratio2 ne "."){$ratio_max = max($ratio1,$ratio2);	}
	if(!exists $disease_ratio{$disease}){print "$disease\t$site_info\n";}
	if($ratio_max < sqrt($disease_ratio{$disease}) && $adar=~/AR/i && $adar!~/AD/i )
	{
		$flag_PM2=1;
	}

	my $flag_BS1=0;
	$ratio_max = 0;
	if($ratio1 ne "." && $ratio2 ne "."){$ratio_max = max($ratio1,$ratio2); }
	my $ratio_min;
	if( $ratio1 ne "." && $ratio2 ne "." )
	{
		$ratio_min = min($ratio1,$ratio2);
		if($ratio_min > ( 1-sqrt(1 - $disease_ratio{$disease}) ) && $adar=~/AD/i )
		{
		#	push @acmg,"BS1";
			$flag_BS1=1 if($ratio_max >= 0.001);
		}

		if($ratio_min > sqrt($disease_ratio{$disease}) && $adar=~/AR/i && $adar!~/AD/i )
		{
		#	push @acmg,"BS1";
			$flag_BS1=1 if($ratio_max >= 0.001);
		}
	}
	if($flag_PM2==1 && $flag_BS1==0){push @acmg, "PM2"; }
	if($flag_PM2==0 && $flag_BS1==1){push @acmg, "BS1"; }

	return @acmg;
}


sub pep_trans{ #wph add 190122
	my $p=shift;
	if($p=~/p\.([A-Y]+)(\d+)([A-Y]+)(.*)/g){
		my $locaa=$2;
		my $firstaa =join"",(map{$protein{$_}//'-'}(split//,$1));
		my $secondaa=join"",(map{$protein{$_}//'-'}(split//,$3));
		my $lastaa=$4;
		$p="p\.${firstaa}${locaa}${secondaa}${lastaa}";
	}elsif($p=~/p\.([A-Y]+)(\d+)fs/g){
		my $locaa=$2;
		my $firstaa =join"",(map{$protein{$_}//'-'}(split//,$1));
		$p="p\.${firstaa}${locaa}fs";
	}elsif($p eq '-'){
		return $p;
	}else{  
		return $p;       
	}
	return $p;
}

sub judge_lastest_sitedb{
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
