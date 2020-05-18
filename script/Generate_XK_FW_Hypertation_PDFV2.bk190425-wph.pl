#!/usr/bin/perl
use strict;
use warnings;
use feature qw/say/;
use File::Find;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw/max min sum maxstr minstr shuffle/;
use Getopt::Long;

my ($final_2_path,$step,$out,$public,$case_detection,$hyper_sample_info_path,$gene_desease_path,$control_detection_path,$disease_file,$latest_gxybase_path,$latest_gxybase_path_fubiao,$latest_gxyreportbase_path,$latest_deletion_path,$hypers_cfg_addsex_path,$NOVO_1256_AF,$HYper_new_AF,$latest_gxybase_old_path);
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
$latest_gxybase_path ||="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY";  #高血压位点库 #zhouweichange 20170629
$latest_gxybase_path_fubiao ||= "/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY_fubiao/"; 
$latest_gxyreportbase_path ||="/PUBLIC/pipline/database/knowledge_database/sitedatabase/GXY_report_base/"; #高血压出报告位点库 #sunxiaoqing 180511
$latest_gxybase_old_path ||= "/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY//201711221646/XK_gxybase201711221646.right.txt";  #高血压位点库 #sunxqchange 20171130
$latest_deletion_path ||= "/PUBLIC/pipline/database/knowledge_database/sitedatabase/filter_site"; # 待删除库
$hypers_cfg_addsex_path ||="$public/gxy_all_cfg_addsex"; #实时更新的高血压config文件，和对应的gxy_detection_ratio对应的人数一致，如果不一致，请检查mergvcf是否执行成功
$NOVO_1256_AF ||="$public/AF/1256WES_Novogene2016_AF.txt"; #novo control 1000WES等位基因频率
$HYper_new_AF ||="$public/AF/final_head_site_all.frq.out.base.right"; #高血压累积case等位基因频率，目前未实时更新，需要今后完善

#############################主脚本区域###########################################
my (%case_ditacion_ratio, %control_ditaction_ratio, %doctor_patient_info, %genelist, %gene_disease, %disease, %hash_base, %hash_base2 ,%hash_base_fubiao, %hash_base_fubiao2 ,%database_site ,%hash_base_2,%hash_del,%right_old, %sampleid_sex, %novo_1256_af, %hyper_last_af);

our %protein=(
	'A'=>'Ala','C'=>'Cys','D'=>'Asp','E'=>'Glu','F'=>'Phe','G'=>'Gly','H'=>'His','I'=>'Ile','K'=>'Lys','L'=>'Leu',
	'M'=>'Met','N'=>'Asn','P'=>'Pro','Q'=>'Gln','R'=>'Arg','S'=>'Ser','T'=>'Thr','V'=>'Val','W'=>'Trp','Y'=>'Tyr',
	'X'=>'*',
);  #wph add 190122

my @titel=(
"Sample",
"CHROM\tStart\tEnd\tREF\tALT\tFunc\tGene\tGeneDetail\tExonicFunc\tAAChange\tNo_annovar_trans\ttransvar_Aachange\tpdf_NM_ID\tpdf_CDS\tpdf_PEP\tgenomicSuperDups\trepeat\tavsnp150\tRare_SNP\tgnomAD_exome_EAS\t1000g2015aug_eas\texac2015_03newEAS\tCombineWES\tphyloP20way_mammalian\tphyloP20way_mammalian_rankscore\tphastCons20way_mammalian\tphastCons20way_mammalian_rankscore\tInterpro_domain\tada_score\trf_score\tHET/HOM\tdpsi_max_tissue\tdpsi_zscore\tProtein\tSIFT;Polyphen2_HVAR;Polyphen2_HDIV;M-CAP;Revel\tSIFT\tPolyphen2_HDIV\tPolyphen2_HVAR\tMCAP\tMutationTaster/CADD_phred\tcytoBand\tGENE_strand\tWarning\tTYPE\tsnp/indel\tDP\tALT_rate\tcds_stop_dis\tbed_edge_dis\tgnomad_depth\tuniprot_domain\tgene_disease\tsite_hotspot\tOMIM\tgene_HGMD_pathogen_contex\tgene_Clivar_pathogen_contex\tHGMD_pathogen_PVS1\tClivar_pathogen_PVS1\tHGMD_pathogen_15bp\tClivar_pathogen_15bp\tHGMD_True_Flase\tClinvar_True_Flase\thgmd_true_genome_contex\thgmd_true_pep_contex\tclivar_true_genome_contex\tclivar_true_pep_contex\tVCF_INFO\tVCF_FORMAT\tVCF_sampleID\tchr_hg38\ts_hg38\te_hg38\tref_hg38\talt_hg38",
#wph changed hgmd,clinvar 190329
#"CHROM\tStart\tEnd\tREF\tALT\tFunc\tGene\tGeneDetail\tExonicFunc\tCDS\tPEP\tAAChange\tcytoBand\tavsnp150\tRare_SNP\tgnomAD_genome_ALL\tgnomAD_exome_ALL\tgnomAD_exome_EAS\tgnomAD_exome_SAS\t1000g2015aug_all\t1000g2015_Chinese\t1000g2015_exac03newALL\t1000g2015_exac03newEAS\t1000g2015_CombineWGS\t1000g2016_CombineWES\tDP\tDP4\tDP4/DP\tALT_RATE\tMQ\tHET/HOM\tWarning\tStrand prefrence\tProtein\tSIFT;Polyphen2_HVAR;Polyphen2_HDIV;M-CAP;Revel\tInterpro_domain\tphyloP20way_mammalian\tphyloP20way_mammalian_rankscore\tphastCons20way_mammalian\tphastCons20way_mammalian_rankscore\tada_score\trf_score\tdpsi_zscore\tTYPE\tsnp/indel\tHGMD;COSMIC;Clinvar\tAA_type\tcosmic81\tclinvar_type\tchr_DB\ts_DB\te_DB\tGene_DB\tHGMD_ID_DB\tAachange_DB\tDisease_name\tclass\tChina/English/inheritance\tpaper\tPubMed\tPhenotypes\tpubmed_id\tMim_Number\tGene_Symbols\tomim_id\tsymbol\tname\talias_symbol\tccds_id\tchr_hg38\ts_hg38\te_hg38\tref_hg38\talt_hg38",
"Disease_system",
"ACMG条目\t致病性结论\t疾病描述\t基因描述\t位点描述\t参考文献\t疾病名称\t遗传模式\t基因名称\tCDS\tPEP\t临床指导建议\t临床表现\tpdf_NMID",
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
&sub_hash_base;
&sub_hash_fubbase;
&sub_hash_oldbase;
&sub_sampleid_sex;
&sub_case_ditacion_ratio;
&sub_control_ditaction_ratio;
&sub_novo_1256_af;
&sub_hyper_last_af;


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
        ($linen,my @site_dase_info)=&add_sitebase_info($linen);       #添加解读库信息
        my @Disease_Namein=(&get_jiyinmiaoshu($cut[6]))[7,1];    #更新筛选的疾病系统

        my $cont_ditacion_term=$control_ditaction_ratio{join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-',1256); #NOVO1256样本检出率
        my $case_ditacion_term=$case_ditacion_ratio{join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-','-');     #高血压样本检出率（实时更新）
#       say "订单${sample_id}中${key}位点未在高血压位点检出率文件（路径为：${control_detection_path}）中存在，请核查${sample_id}样本是否执行过mergevcf.sh程序。如未执行，请执行此程序；如有执行，请核查高血压config文件（路径：${hypers_cfg_addsex_path}）最新的文件中是否存在此订单，而高血压检出率文件中未找到此订单。如有此种情况，请核查mergevcf.sh的正确，并在高血压config文件去掉此订单，并重新运行mergevcf.sh程序" if ($case_ditacion_term eq "-\t-\t-\t-");
        
        my @all_sample_case=split/,/,((split/\t/,$case_ditacion_term)[2]);
        my @all_sample  =map{$doctor_patient_info{$_}->[0] // '-'}@all_sample_case;
        my @all_dotcor  =map{$doctor_patient_info{$_}->[1] // '-'}@all_sample_case;
        my @all_hospital=map{$doctor_patient_info{$_}->[2] // '-'}@all_sample_case;

        my @all_sample_gender=map{$sampleid_sex{$_} // '-'}@all_sample_case;
        my $patient_doctor_hos=(@{$doctor_patient_info{$sample_id}})  ? $doctor_patient_info{$sample_id} : [("-")x3]; #医生患者信息
        my $patient_sex_inform=$sampleid_sex{$sample_id} // "-";  #性别信息
        my $case_control_rates=
           ((split/\t/,$case_ditacion_term)[0]=~/[0-9]+\.[0-9]*/i and (split/\t/,$case_ditacion_term)[0]>0 and
            (split/\t/,$cont_ditacion_term)[0]=~/[0-9]+\.[0-9]*/i and (split/\t/,$cont_ditacion_term)[0]>0 )
           ?(split/\t/,$case_ditacion_term)[0]/(split/\t/,$cont_ditacion_term)[0] : '-' ;              #高血压样本检出率 和 NOVO1256样本检出率 比例

        my $hyper_last_af_case=$hyper_last_af{join"\t",@cut[0,1,3,4]} // join"\t",(("-")x4); # NOVO1256样本AF
        my $novo_1256_af_contr=$novo_1256_af{join"\t",@cut[0,1,3,4]} // join"\t",(("-")x4);  # 高血压样本AF（未实时更新，需要增加）
        #say "订单${sample_id}中${key}位点未在高血压阳性样本AF文件（$HYper_new_AF），文件未实时更新。目前最新更新市建委20170416" if ($novo_1256_af_contr eq "-\t-\t-\t-");

        my @final_out= (
            $sample_id,$linen,$diseaseasystem,@site_dase_info,@Disease_Namein,
            @$patient_doctor_hos,$patient_sex_inform,(join",",@all_sample),(join",",@all_dotcor),(join",",@all_hospital),(join",",@all_sample_gender),
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
        for(@cut[19..22]){($_ eq "." || $_ <= 0.01) && ++$ratio1000};
	#if($cut[16] eq "." || $cut[16] < 0.01)	{++$ratio1000;}
	#if($cut[17] eq "." || $cut[17] < 0.01)  {++$ratio1000;}

        my $sample_detected_num_1256= (split/\t/,$cont_ditacion_term)[1];
        my $sample_detected_num_case= (split/\t/,$case_ditacion_term)[1];
        my @protoin_predict=split/;/,$cut[34];
	my $protoin_predict_num=0;
        my ($SIFT_Polyphen2,$CADD)=(@protoin_predict>=4) ? ( (join"",@protoin_predict[0..3]),$protoin_predict[-1] ) : ($cut[34],$cut[34]); ##changed by sunxq 20171122
        
        ##############过滤###############################
        my @filter_type=();
	my $hgmd_flag=0;
	my $Disease_Namein_flag=0;
	($cut[45]<4)  && (push @filter_type, 'DP<4');
        #($dp4_ratio<0.1) && (push @filter_type, 'DP4_ratio<0.1');
        ($ratio1000<4) && (push @filter_type, 'AF1000_4>0.01');
        ($Disease_Namein[1] eq '-') && (push @filter_type, "No_disease")&&($Disease_Namein_flag=1);
        ($patient_sex_inform=~/男/ and $site_dase_info[6]=~/妊娠高血压综合/) && (push @filter_type, 'Man_pregnancy_hypertension') ;
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
			if($flag_dm=~/DM|pathogenic|significance/i )
			{
				$FLAG_DM_62=2 if($i==62);
				$FLAG_DM_64=2 if($i==64);
			}
		}
	}
	
	if( ($FLAG_DM_T62==1 && $FLAG_DM_T64==1 && $FLAG_DM_62==1 && $FLAG_DM_64==1) || ($FLAG_DM_T62==1 && $FLAG_DM_T64==0 && $FLAG_DM_62==1 && $FLAG_DM_64==0) || ($FLAG_DM_T62==0 && $FLAG_DM_T64==1 && $FLAG_DM_62==0 && $FLAG_DM_64==1) )
	{
		push @filter_type, "non_pathogenic_DB" if($mu_type ne "synonymous SNV");
	}else{
		$hgmd_flag=1 if($mu_type ne "synonymous SNV");
	}

	if( ($FLAG_DM_T62==1 && $FLAG_DM_62==2) || ($FLAG_DM_T64==1 && $FLAG_DM_64==2) )
	{
		$hgmd_flag=1 if($mu_type eq "synonymous SNV");
	}else{
		push @filter_type, "synonymous_SNV" if($mu_type eq "synonymous SNV");
	}
	
	#$hgmd_flag=1 if( $cut[43]=~/true/i && $cut[43]=~/DM/i );  ### panqi 190124
	
	my $flag=(@filter_type) ? (join",",@filter_type) : '?';
	if (($cut[52] ne "-" and $flag!~/Man_pregnancy_hypertension/i and $cut[52]!~/negative/ ) or ($hgmd_flag==1 && $ratio1000>2 && $Disease_Namein_flag==0 && @filter_type<1)  or @filter_type<1){  ### panqi 181017 add ratio1000 判断
		++$sample_number{$sample_id};
		if(exists $hash_base{$key2})
		{
			$flag="done";
			$flag.="_NM" if($hash_base2{$key2} eq "-");
		}
		else{	$flag="todo"; }
            	#$flag=(exists $hash_base{$key2}) ? "done" :"todo";
		
		if(exists $hash_base_fubiao{$key2} && !exists $hash_base{$key2})
		{
			$flag="done_fubiao";
			$flag.="_NM" if($hash_base_fubiao2{$key2} eq "-");
		}
		#$flag="done_fubiao" if(exists $hash_base_fubiao{$key2} && !exists $hash_base{$key2});
		
		$flag="todo_oldDB" if(exists $right_old{join"\t",@cut[0,1,3,4]}  && !exists $hash_base{$key2} && !exists $hash_base_fubiao{$key2});  #######sunxq 20171205
		if($hgmd_flag==1){
			if(($cut[52] ne "-" and $flag!~/Man_pregnancy_hypertension/i) or @filter_type<1){
			}else{
				$flag.="_DM";
			}
		}
            	say RIGHT "$flag\t$all_output";
            	#if ($flag eq "base_to_read"){
                #	say READ  join"\t",("-\t-\t-",(split/\t/,$case_ditacion_term)[2],$_,$diseaseasystem,@site_dase_info) ;
		#}
        }else{
                if(exists $hash_base{$key2}){$flag.=";exist_in_readbase";}
		if(exists $hash_base_fubiao{$key2} && !exists $hash_base{$key2}){$flag.=";exist_in_FubiaoReadbase";}
		if($Disease_Namein_flag==0 || $cut[43]=~/true_positive/i){
            		say WRONG "$flag\t$all_output";                    #### panqi 180406 ???
		}
        }
    }
    close IN;
}
foreach my $idsamle (keys%sample_number){
    ($sample_number{$idsamle}<=2) && (say RIGHT join"\t",('Negative_sample',"$idsamle",('-')x103));  #### panqi 180406 ???
}
say "step 1执行完成！！";
close FI;
close ALL;
close RIGHT;
close WRONG;
}

sub sub_step_2{
#&sub_hash_base;
#&sub_hash_deletion;
my %sampledata=&read_samdatabase;
my %coresitedb;
my %sample_flag;

open RIGHT,   ">>$out/GXY_sample_rigth_info.txt"   or die "Cannot find the File :$out/sample_rigth_info.txt";
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


open RIGHT,   "$out/GXY_sample_rigth_info.txt"   or die "Cannot find the File :$out/sample_rigth_info.txt";
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
	say FINALPDF join"\t",@cutnew[1..70,76,71..75,55,77..89,90,-2];#wph add 解读人20181207
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
close CORESITE;
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
########################## 高血压位点库文件 ###################################################
sub sub_hash_base{
    my $latest_gxybase_prefix=".right.txt";
    my $latest_gxybase=&get_latest_data($latest_gxybase_path,$latest_gxybase_prefix);
say '*'x25;
say "正在读取最新高血压位点库文件： 
$latest_gxybase";
    open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
        %hash_base=map{chomp; (join"\t",(split/\t/,$_)[3,4,6,7] ) => (join"\t",(split/\t/,$_)[78..90]) ;}<DATA>;
    close DATA;
        open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
        %database_site=map{chomp; (join"\t",(split/\t/,$_)[3,4,6,7]) => (join"\t",(split/\t/,$_)[11..14]) ;}<DATA>;
    close DATA;
    open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
	%hash_base2=map{chomp; (join"\t",(split/\t/,$_)[3,4,6,7] ) => (join"\t",(split/\t/,$_)[91]) ;}<DATA>;
    close DATA;
    open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
	%hash_base_2=map{chomp; (join"\t",(split/\t/,$_)[2,3,4,6,7]) => 1}<DATA>;
    close DATA;
say '读取完成！！';
}

##################高血压附表位点库########3
sub sub_hash_fubbase{
	my $latest_gxybase_prefix=".right.txt";
	my $latest_gxybase=&get_latest_data($latest_gxybase_path_fubiao,$latest_gxybase_prefix);
say '*'x25;
say "正在读取最新高血压附表位点库文件：
$latest_gxybase";
	open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
		%hash_base_fubiao=map{chomp; (join"\t",(split/\t/,$_)[3,4,6,7] ) => (join"\t",(split/\t/,$_)[78..90]) ;}<DATA>;
	close DATA;
	open DATA,$latest_gxybase or die "can't open the file:$latest_gxybase";
		%hash_base_fubiao2=map{chomp; (join"\t",(split/\t/,$_)[3,4,6,7] ) => (join"\t",(split/\t/,$_)[91]) ;}<DATA>;
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
                $hash{$keysam}{info}=join("\t",@tmp[2..$n]);
                $hash{$keysam}{key}=$tmp[0]; 
        }
        close  DATA;
        return %hash; 
say '读取完成！！';
}
##########################################################
sub sub_hash_oldbase{
say '*'x25;
say "正在读旧版取高血压位点库文件：
$latest_gxybase_old_path";
	open DATA,$latest_gxybase_old_path or die "can't open the file:$latest_gxybase_old_path";
		%right_old=map{chomp; (join"\t",(split/\t/,$_)[4,5,7,8]) => (join"\t",(split/\t/,$_)[67..79]) ;}<DATA>;
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

sub add_sitebase_info{
    my $site_info=shift;
    next if /^CHROM/;
    my @cut=split/\t/,$site_info;
            
    my @add_info=(("-")x13);
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
    
    if (exists $hash_base{$key}){
        @add_info=split/\t/,$hash_base{$key};
	@add_info[6,2,7,11,3,12]= (&get_jiyinmiaoshu($cut[6]))[1,2,3,4,5,6];
	$add_info[13]=$cut[12]; #CDS    PEP
	if($cut[30]=~/het/ and $add_info[7]=~/AR/){$add_info[4]="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，".$add_info[4];}#wph add 190104
    }elsif(exists $hash_base_fubiao{$key}){
	@add_info=split/\t/,$hash_base_fubiao{$key};
	@add_info[6,2,7,11,3,12]= (&get_jiyinmiaoshu($cut[6]))[1,2,3,4,5,6];
	$add_info[13]=$cut[12]; #CDS    PEP
	if($cut[30]=~/het/ and $add_info[7]=~/AR/){$add_info[4]="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，".$add_info[4];}#wph add 190104
    }else{
	my $extra_type="";
        my $mutation_type=
            ($cut[8]=~/stopgain/i)      ? "无义"  :
	    ($cut[8]=~/stoploss/i)      ? "终止缺失"  :
            ($cut[8]=~/splicing/i)      ? "剪切":	
            ($cut[8]=~/nonframeshift/i) ? "非移码":
            ($cut[8]=~/frameshift/i)    ? "移码"  :
            ($cut[8]=~/nonsynonymous/i) ? "错义"  :
            ($cut[8]=~/synonymous/i)    ? "同义"  : "?";
        my $detection_or_not=
            ($cut[19]eq "." and $cut[20]eq ".") ? "均未检出" :'均低于1%';
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
		if($cut[28] > 0.6 and $cut[29] > 0.6 and abs($cut[32]) >=2)
		{
			$splicing_pred.="dbscSNV1.1、SPIDEX软件预测该变异对基因或基因产物有害。";
		}
		elsif($cut[28] > 0.6 and $cut[29] > 0.6)
		{
			$splicing_pred .="dbscSNV1.1软件预测该变异对基因或基因产物有害。";
		}
        elsif(abs($cut[32]) >=2)
        {
            $splicing_pred.="SPIDEX软件预测该变异对基因或基因产物有害。";
        }
        else
        {
            $splicing_pred.="";
        }
	}
	my $wenxianzhichi=
           ($cut[62] eq "-" and $cut[63] eq "-" and $cut[64] eq "-" and $cut[65] eq "-" ) ? "由于暂无相关研究，" : "由于研究较少，";
	#($cut[58]=~/\S/ and (length $cut[58]) >10 and $cut[59]=~/http/) ? "由于研究较少，" : "由于暂无相关研究，";
        
	my $muty="";#wph changed 190104 "$cut[14]${mutation_type}";
	if($cut[5]=~/splicing/){
		$muty="位于剪切区";
	}else{
		$muty="会导致$add_info[10]${mutation_type}变异";
	}
	$muty=~s/^-//;
	if(!exists $right_old{$key}){ 
        	@add_info[0,1]=($cut[43]=~/true/i) ? ("待解读","待解读") : ("自动","临床意义未明") ; #ACMG条目,#致病性结论
        	@add_info[8,6,2,7,11,3,12]= (&get_jiyinmiaoshu($cut[6]))[0..6];  #致病基因 疾病名称 疾病描述 遗传模式 临床指导建议 基因描述    临床表现 疾病大类
        	@add_info[9,10,13]=@cut[13,14,12]; #CDS    PEP
		($add_info[10],$extra_type)=&pep_trans($cut[14]); #wph add 190122
		if($mutation_type eq "?" or $cut[8] =~/unknown/i){ #wph add 190418
			$mutation_type=$extra_type;
			$cut[8] =
			($extra_type eq "错义")					? "nonsynonymous SNV"	    :
			($extra_type eq "同义")					? "synonymous SNV"	    :
			($extra_type eq "无义")                                 ? "stopgain"                :
			($extra_type eq "终止缺失")                             ? "stoploss"                :
($extra_type eq "剪切" && $cut[13] !~/del|->|>-/i && $cut[5] ne "intronic")	? "splicing SNV"	    :
			($extra_type eq "剪切" && $cut[13] =~/(del|->|>-)/i)	? "splicing INDEL"	    :
			($extra_type eq "非移码" && $cut[13] =~/del/) 		? "nonframeshift deletion"  :
			($extra_type eq "非移码" && $cut[13] =~/dup/) 		? "nonframeshift insertion" :
			($extra_type eq "移码" && $cut[13] =~/del/)		? "frameshift deletion"	    :
			($extra_type eq "移码" && $cut[13] =~/dup/)		? "frameshift insertion"    : "SNV";
		}
		if($cut[5]=~/splicing/ or $extra_type eq "剪切"){
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
		my $arhet="";#wph changed 190104
		if($cut[30]=~/het/ and $add_info[7]=~/AR/){$arhet="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，";}
        	$add_info[4]="本筛查检测出$cut[6]基因的变异，变异位点为$cut[13]，查询ClinVar等公共数据库显示，$cut[13]变异${muty}。$ratio_detection"."${software_pred}${splicing_pred}${wenxianzhichi}行ACMG标准，判定该变异的致病性为临床意义未明。${arhet}具体情况请结合临床相关资料综合判断。";#位点描述
#该突变在千人基因组和ExAC数据库中的基因频率${detection_or_not}"."${software_pred}。${wenxianzhichi}行ACMG标准，判定该突变的致病性为临床意义未明。具体情况请结合临床相关资料综合判断。";#位点描述
        	$add_info[5]="-"; #参考文献
	}

	if(exists $right_old{$key}){                                       ### panqi 20181025
		@add_info=split/\t/,$right_old{$key};
		my $site_info_tmp=$add_info[4];
		@add_info[8,6,2,7,11,3,12]= (&get_jiyinmiaoshu($cut[6]))[0,1,2,3,4,5,6];
		$add_info[13]=$cut[12]; #CDS    PEP
		my $arhet="";#wph changed 190104
		if($cut[30]=~/het/ and $add_info[7]=~/AR/){$arhet="由于该疾病为隐性遗传病，本次检出杂合变异，一般情况下不会出现疾病的相应症状或者症状轻微，";}
		$add_info[4] = "本筛查检测出$cut[6]基因的变异，变异位点为$cut[13]，查询ClinVar等公共数据库显示，$cut[13]变异${muty}。${ratio_detection}"."${software_pred}${splicing_pred}${wenxianzhichi}行ACMG标准，判定该变异的致病性为临床意义未明。${arhet}具体情况请结合临床相关资料综合判断。"."<-->$site_info_tmp";        #位点描述
		#$add_info[4] .= " \( $site_info_tmp \) ";                 ### panqi 失效
	}
    }
	$site_info=join("\t",@cut);
    return $site_info,@add_info;
}
close IN;

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
	}elsif($p eq '-'){
		$type="剪切";
		return $p,$type;
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

