package StandardClass;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(  
        Standard_confirm,acmg_auto_judge
);



################################## sub Meet standards
sub Standard_confirm{
	my ($hash,$genehash,$sampleID)=@_;
	my($flag_exists,$flag_DP,$flag_Rate,$flag_protein,$flag_dm,$flag_benign,$flag_hot,$flag_Negative,$mu_type,$flag_snp,$flag_synpath)=(0,0,0,0,0,0,0,0,0,0,0);
	
	#####判断基因是否在检测列表中
	my @allg=split /;/,$$hash{'Gene'};
	foreach my $k(@allg){
		$flag_exists=1 if(exists $$genehash{$sampleID}{$k} && $$genehash{$sampleID}{$k}==1);  #在基因列表中，标签值为1
	}

	#####判断是snp还是indel
	$flag_snp=1 if($$hash{"snp/indel"} =~/snp/i);
	#####判断深度

	$flag_DP=1 if($$hash{"DP"} >=10);
	
	#####判断变异类型
	$mu_type=1 if($$hash{"ExonicFunc"} eq "synonymous SNV" || ($$hash{"Func"}!~/splicing/ && $$hash{"Func"}!~/exonic/) );    #"外显子/剪切区变异，同义突变为1，非外显子剪切区变异为1，外显子剪切区变异为0（错义，stopgain，stoploss，splicing，unknown）"
	
	####判断频率
	my $Ratenum=0;
	my @arrRatio=("gnomAD_exome_EAS","exac2015_03newEAS","CombineWES"); # gnomAD_exome_EAS	1000g2015aug_eas(这个不看)	exac2015_03newEAS	CombineWES

	for(@arrRatio){
		next if($$hash{$_} eq "." || $$hash{$_} <=0.001);
		$Ratenum++;
	}
	$flag_Rate=1 if($Ratenum==0);  #满足频率要求，标签值为1

	#####判断蛋白软件有害性
	my @protoin_predict=split/;/,$$hash{"SIFT;Polyphen2_HVAR;Polyphen2_HDIV;M-CAP;Revel"};
	my $protoin_predict_num=0;
	my $flag_BP7=0;
	my $splicing_syn_tag=0;
	for(@protoin_predict){
		if($_=~/tolerated|benign/i){$protoin_predict_num++;}
		if($_=~/deleterious|damaging/i){$protoin_predict_num=-1;} #wph add 190115
	}
	

	if( $$hash{"ExonicFunc"} =~/^synonymous/ || $$hash{"Func"}=~/^splicing/ )
	{
		$flag_BP7++ if($$hash{"ada_score"} eq "." );
		$flag_BP7++ if($$hash{"rf_score"} eq "."  );
		$flag_BP7++ if($$hash{"dpsi_zscore"} eq "." );
		$splicing_syn_tag=1;
	}
	
	if($splicing_syn_tag==1){
		$flag_protein=1 if($flag_BP7<3); #对于剪切区域和同义突变的位点，如果这两个软件都没有预测出来就放error
	}else{
		$flag_protein=1 if($protoin_predict_num<3) #四个良性或者3个良性1个-算一类标签，则蛋白软件不通过，其余全通过。
	}

	###########良性、有害位点判断
	my @truearry=("hgmd_true_genome_contex","clivar_true_genome_contex");
	my $pathNum=0;
	my $flag_clinvar=0;
	my $beignNum=0;
	foreach  (@truearry) {
		my $datab=$_;
		next if($$hash{$datab} eq "-");
		my @flag_dm_array=split /\|/,$$hash{$_};
		for(@flag_dm_array){
			$flag_dm_type=(split /:/,$_)[0];
			next if($flag_dm_type!~/true_type/);
			my $flag_dm=(split /:/,$_)[3];
			if($flag_dm!~/DM|pathogenic|significance/i && $flag_dm ne "not provided" && $flag_dm ne "other" ){
				$beignNum++;
				$flag_clinvar=1 if ($datab =~ /clivar/i);
			}
			if($flag_dm=~/DM|pathogenic/i  && $flag_dm!~/conflict/i ){
				$pathNum++;
			}
		}
	}
	$flag_dm=1 if($pathNum>0); #至少有一个为type致病为致病
	$flag_benign=1 if ($beignNum>1) ; #两个都为 type 良性的点为良性点
	##############热点、删除的点
	if($$hash{"site_hotspot"} ne '-'){
		if($$hash{"site_hotspot"}=~/negative/i){
			$flag_Negative=1; 
		}else{
			$flag_hot=1;
		}
	}
	#######同义突变基因
	if($$hash{"gene_Clivar_pathogen_contex"} =~/synoonymous:(\d+)/){
		my $num=$1;
		$flag_synpath=1 if($num>0);
	}
	######
	my $errorTags=&getTag($flag_exists,$flag_DP,$flag_Rate,$flag_protein,$flag_dm,$flag_benign,$flag_hot,$flag_Negative,$$hash{"Func"},$$hash{"ExonicFunc"},$flag_clinvar);
	#######
	return $flag_exists,$flag_DP,$flag_Rate,$flag_protein,$flag_dm,$flag_benign,$flag_hot,$flag_Negative,$mu_type,$flag_snp,$flag_synpath,$errorTags;
}

##########
sub getTag{
	my ($flag_exists,$flag_DP,$flag_Rate,$flag_protein,$flag_dm,$flag_benign,$flag_hot,$flag_Negative,$func,$ExonicFunc,$flag_clinvar)=@_;
	my @tags=();
	push @tags,"DP<10" if($flag_DP==0);
	push @tags,"AF" if($flag_Rate==0);
	push @tags,"protion_predict_filter" if($flag_protein==0);
	push @tags,"Del_site_hotspot" if($flag_Negative==1);
	push @tags,"intronic" if($func=~/intronic/);
	push @tags,"UTR" if($func=~/UTR/ );
	push @tags,"synonymous_SNV" if( $ExonicFunc eq 'synonymous SNV' );
	push @tags,"Clinvar_Benign" if($flag_clinvar==1);
	my $jointags=join(";",@tags);
	return $jointags;
}


sub acmg_auto_judge{                    #wph add 190312       
	my $hash  = shift;
	my $sample  = shift;
	my $site_info  = shift;
	my %clivarPath = shift;
	my $flag_Rate = shift;
	my $control_ditaction_ratio = shift;
	my $case_ditacion_ratio = shift;

	my @cut=split/\t/,$site_info;
	my @acmg=();
	my ($pep1,$pep2,$pep3)=("","","");
	my @site;

########### 蛋白软件预测
	if($$hash{'Protein'} == 1 && $$hash{"ExonicFunc"}=~/^nonsynonymous/ ){
                push @acmg,"PP3";  
	}
	
	my $flag_BP7=0;
	if( $$hash{"ExonicFunc"} =~/^synonymous/ || $$hash{"Func"}=~/^splicing/ )
	{
		
		$flag_BP7++ if($$hash{"ada_score"} ne "." && $$hash{"ada_score"}!~/SNV/ && $$hash{"ada_score"}!~/exon/ && $$hash{"ada_score"} > 0.6);
		$flag_BP7++ if($$hash{"rf_score"} ne "." && $$hash{"rf_score"}!~/SNV/ && $$hash{"rf_score"}!~/exon/ && $$hash{"rf_score"} > 0.6);
		$flag_BP7++ if($$hash{"dpsi_zscore"} ne "." && abs($$hash{"dpsi_zscore"} ) >= 2 );
		if($flag_BP7==3 )
		{
			push @acmg,"PP3";  
		}
	}
###############数据库1
	my $flag_PS1=0;
	my $flag_truetype=0;
	my %star=('one'=>1,'two'=>2,'three'=>3,'four'=>4,'none'=>'0');
	my $clinvar_alle='';
	my ($truehgmd)= $$hash{"HGMD_True_Flase"}=~/genome:([^\|]+)/;
	my ($truecli)=$$hash{"Clinvar_True_Flase"}=~/genome:([^\|]+)/;
	my ($truecliPep)=$$hash{"Clinvar_True_Flase"}=~/pep:([^\|]+)/;
	if( $$hash{"ExonicFunc"}=~/^nonsynonymous/ && $truehgmd!~/true_type/i  && $truecli!~/true_type/i && $truecliPep=~/true_type/)
	{
		@site=split /\|/,$$hash{"clivar_true_pep_contex"};
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
	if($flag_PS1==1 && $flag_truetype==1 && $flag_Rate==1){ push @acmg,"PS1"; print "PS1\t$site_info\n";}
#############数据库2
	my $flag_PM5=0;
	$flag_truetype=0;
	my $flag_Ter=0;
	
	if( ( $$hash{"ExonicFunc"}=~/^nonsynonymous/) && $$hash{"clivar_true_pep_contex"} ne "-" ) 
	{
		@site=split /\|/,$$hash{"clivar_true_pep_contex"};
		for(@site)
		{
			$flag_truetype=1 if( (split /:/,$_)[0] =~/true_type/i );
			next if( (split /:/,$_)[0] !~ /true_site/i );
			$pep2=(split /:/,$_)[3];
			if( $pep2=~/^Pathogenic/ || $pep2=~/^Likely pathogenic/ )
                        {
                                $flag_PM5=1;
                                if( (split /:/,$_)[2]=~m/p\.(\S+)(\d+)(\S+)/ || (split /:/,$_)[2]=~m/p\.(\S)(\d+)(\*)/ )
                                {
                                        $flag_Ter=1 if($3=~"Ter" || $3 eq "*" || $3 eq "X");
                                }
                        }
		}
	}
	if($flag_PM5==1 && $flag_truetype==0 && $flag_Ter==0 && $flag_Rate==1){ push @acmg,"PM5"; }

	
###############################频率
	my $flag_PM2_1=0;
	my $flag_PM2_2=0;
	my $db_dep_flag20=0;
	my @db_dep=split /\/\//, $$hash{"gnomad_depth/exac_depth"};
	my ($db_dep_1, $db_dep_2, $db_dep_12, $db_dep_22);
	my ($db_dep_flag1,$db_dep_flag2) = (split /;/, $db_dep[2])[0,1];
	if($db_dep_flag1 eq "Y" )
	{
		$db_dep_1 = (split /\|/, $db_dep[0])[1];
		my @aa=split /,/,$db_dep_1;
		my $numaa=@aa;
		my $depaa=0;
		foreach my $dep (@aa) {
			$depaa+=$dep;
		}
		$db_dep_12= $depaa/$numaa;
	}
	if($db_dep_flag2 eq "Y" )
	{
		$db_dep_2 = (split /\|/, $db_dep[1])[1];
		my @aa=split /,/,$db_dep_2;
		my $numaa=@aa;
		my $depaa=0;
		foreach my $dep (@aa) {
			$depaa+=$dep;
		}
		$db_dep_22= $depaa/$numaa;
	}
	if($db_dep_flag1 eq "Y" && $db_dep_flag2 eq "Y" )
	{
		$db_dep_flag20=1 if( $db_dep_12>=20 && $db_dep_22>=20 );
	} #此处要求两个数据库深度都大于20
	
	if( ($$hash{"gnomAD_exome_EAS"} eq "." || $$hash{"gnomAD_exome_EAS"}  eq "0") && ( $$hash{"exac2015_03newEAS"} eq "." || $$hash{"exac2015_03newEAS"}  eq "0") && $db_dep_flag20==1 )
	{
		$flag_PM2_1=1;
	}
	
	my $cont_ditacion_term=$$control_ditaction_ratio{ join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-'); #NOVO1256样本检出率
	my $case_ditacion_term=$$case_ditacion_ratio{     join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-');
	if( ( ( (split/\t/,$case_ditacion_term)[1]==1 && ((split/\t/,$case_ditacion_term)[2] =~/$sample$/ || ($case_ditacion_term)[2] =~/$sample\|/) ) || (split/\t/,$case_ditacion_term)[1]==0 )  && (split/\t/,$cont_ditacion_term)[1]==0 )
	{
		$flag_PM2_2=1 if( $$hash{"snp/indel"} eq "snp");
	}
	if($flag_PM2_1==1 || $flag_PM2_2==1){ push @acmg,"PM2";  }

	####################
	
	return $db_dep_flag20,\@acmg;
}

