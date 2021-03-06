#use strict;
#use warnings;
package PathogenicityClass;
use Exporter;
use List::Util qw/max min sum maxstr minstr shuffle/;
our @ISA = qw(Exporter);
our @EXPORT = qw(  
        PathogenDatabase acmg_auto_judge PathogenicityJudge
);

my %star=('one'=>1,'two'=>2,'three'=>3,'four'=>4,'none'=>'0');

sub PathogenDatabase{
	my ($ClivarDatabase,$snpindel,$gene,$trueHGMD,$trueTexHGMD,$trueClivar,$trueTexClivar,$trueTexClivarP,$omimhash,$hgnchash)=@_;
	my ($true1,$true2)=('-','-');
	my ($truehgmd)=$trueHGMD=~/genome:([^\|]+)/;
	my ($truecli)=$trueClivar=~/genome:([^\|]+)/;
	my ($truecliPep)=$trueClivar=~/pep:([^\|]+)/;
	my $judgeResult="null";
	my $geneOmimR="-";
	#print "$gene\t$snpindel\t";
	#print "$trueTexHGMD\n";
	#print "$trueTexClivar\n";
	if($truehgmd!~/true_type/i  && $truecli!~/true_type/i){
		#print "$geneOmimR\t0\t";
		if($snpindel eq 'snp' && $truecliPep=~/true_type/){
			my ($clipath,$alleid)=&checkDM('Clivar',$trueTexClivarP,\%$ClivarDatabase);
			my $starr=$star{$$ClivarDatabase{$alleid}{star}};
			if($clipath eq 'path'){
				if($starr>=2){
					$judgeResult='PS1';
				}
			}
		}
	}elsif($truehgmd=~/true_type/i && $truecli!~/true_type/i){ #print "$geneOmimR\t1\t";
		my ($hgmdpath,$alleid)=&checkDM('HGMD',$trueTexHGMD,\%$ClivarDatabase);
		if($hgmdpath eq 'path'){
			$judgeResult='PP';
		}
	}elsif($truehgmd!~/true_type/i && $truecli=~/true_type/i){
		my ($clipath,$alleid)=&checkDM('Clivar',$trueTexClivar,\%$ClivarDatabase);
		if($$ClivarDatabase{$alleid}{OMIM} eq '-'){
			$geneOmimR='N'	
		}else{	
			$geneOmimR=&geneOmimcheck($alleid,\%$ClivarDatabase,$omimhash,\%$hgnchash);
		}
		my $starr=$star{$$ClivarDatabase{$alleid}{star}};
		# print "$geneOmimR\t2\t";
		if($clipath eq 'path' && $geneOmimR ne 'W'){
			if($starr>=2){
				$judgeResult='可能致病';
			}else{
				$judgeResult='PP';
			}
		}elsif($clipath eq 'benign' && $geneOmimR ne 'W'){
			$judgeResult='BP6';
		}
	}elsif($truehgmd=~/true_type/i && $truecli=~/true_type/i){
		my ($hgmdpath,$alleidh)=&checkDM('HGMD',$trueTexHGMD,\%$ClivarDatabase);
		my ($clipath,$alleidc)=&checkDM('Clivar',$trueTexClivar,\%$ClivarDatabase);
		if($$ClivarDatabase{$alleidc}{OMIM} eq '-'){
			$geneOmimR='N'
		}else{
			$geneOmimR=&geneOmimcheck($alleidc,\%$ClivarDatabase,$omimhash,\%$hgnchash);
		}
		my $starr=$star{$$ClivarDatabase{$alleidc}{star}};
		if( $geneOmimR ne 'W'){  #print "$geneOmimR\t3\t";
			if($hgmdpath eq 'path' && $clipath eq 'path'){
				if($starr>=2){
                                	$judgeResult='可能致病';
                        	}else{
                                	$judgeResult='PP+PP';
                        	}
			}elsif($clipath eq 'benign' && $hgmdpath ne 'path'){
				$judgeResult='BP6';
			}elsif($clipath eq 'path' && $hgmdpath ne 'path'){
				if($starr>=2){
                                	$judgeResult='可能致病';
				}else{
                                	$judgeResult='PP';
				}
			}elsif($clipath eq 'benign' && $hgmdpath eq 'path'){
				$judgeResult='BP6+PP';
			}elsif($clipath eq 'uncertain' && $hgmdpath eq 'path'){
				$judgeResult='PP';
			}elsif($clipath eq 'conflict' && $hgmdpath eq 'path'){
				$judgeResult='PP';
			}
		}else{	# print "$geneOmimR\t4\t";
                	if($hgmdpath eq 'path'){
                        	$judgeResult='PP';
                	}
		}
	}
	#print "$truehgmd\t$truecli\t$trueHGMD\t$trueTexClivar\n";
	$judgeResult.="-$geneOmimR";
	return $judgeResult;
}


sub checkDM{
	my ($database,$baseline,$hash)=@_;
	my $pathR='false';
	my $AlleleID='-';
	my $alleLevel=0;
	if($database eq 'HGMD'){
		$baseline=~s/\|\|/\/\//g;
		my @a=split /\|/,$baseline;
		foreach my $tit(@a){
			my @aaa=split /:/,$tit;
			if($aaa[0]=~/true_type/i && $aaa[3]=~/DM/){
				$pathR='path';
			}
		}
	}elsif($database eq 'Clivar'){
		$baseline=~s/\|\|/\/\//g;
                my @a=split /\|/,$baseline;
		my $flagben=0;
		my $flagpath=0;
		my $flagconflict=0;
		my $flaguncertain=0;
		foreach my $tit(@a){
			my @aaa=split /:/,$tit;
			if($aaa[0]=~/true_type/i){
				my ($pathtit)=(split /,/,$aaa[3])[0];
				if($pathtit=~/^Benign/ || $pathtit=~/^Likely benign/){
					$flagben=1;
				}
				if($pathtit=~/^Pathogenic/ || $pathtit=~/^Likely pathogenic/){
					$flagpath=1;		
				}
				if($pathtit=~/Conflicting interpretations of pathogenicity/i){
					$flagconflict=1;
				}
				if($pathtit=~/Uncertain significance/i){
					$flaguncertain=1;
				}
				my($aid)=(split /_/,$aaa[1])[0];
				if($star{$$hash{$aid}{star}}>=$alleLevel){
					$alleLevel=$star{$$hash{$aid}{star}};
					$AlleleID=$aid;
				}
			}
		}
		if($flagconflict==1){
			$pathR='conflict';
		}elsif($flagpath==1){
			$pathR='path';
		}elsif($flagben==1){
			$pathR='benign';
		}elsif($flaguncertain==1){
			$pathR='uncertain';
		}
	}
	return $pathR,$AlleleID;
}

sub geneOmimcheck{
        my($alleid,$ClivarDatabase,$omimhash,$hgnchash)=@_;
        my @omims=split /;/,$$ClivarDatabase{$alleid}{OMIM};
        my $check='W';
	my ($omimid)=(split /_/,$omimhash)[0];
	foreach my $omim2(@omims){
        	if($omimid eq $omim2){
			$check="R";
        	}
	}
        return $check;
}


sub geneOmimcheck_bk{
	my($gene,$alleid,$ClivarDatabase,$omimhash,$hgnchash)=@_;
	my @genes=split /;/,$gene;
	my @omims=split /;/,$$ClivarDatabase{$alleid}{OMIM};
	my $check='W';
	foreach my $gg(@genes){
		my $g=$$hgnchash{$gg};
		foreach my $omim1(keys %{$$omimhash{$g}}){
			foreach my $omim2(@omims){
				if($omim1 eq $omim2){
					$check="R";
				}
			}
		}
	}
	return $check;
}

sub acmg_auto_judge{                    #wph add 190312       
	my ($site_info,$sample,$cont_ditacion_term,$case_ditacion_term) = @_;
	#my $sample    = shift;
	#my $cont_ditacion_term = shift;
	#my $case_ditacion_term = shift;
	#my $disease   = shift;
	#my $adar      = shift;

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

	my ($pep3_num1,$pep3_num2,$pep3_num3,$pep3_rate,$pep3_flag)=(0,0,0,0.0001,0);
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

	my $flag_PM5=0;
	my $flag_truetype=0;
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
			#$flag_PM5=1 if( $pep2=~/Pathogenic/ || $pep2=~/Likely pathogenic/ );
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

	$pep3_flag=0;
	if( $cut[59] ne "-"  &&   $cut[8]=~/^nonsynonymous/)
	{
		
		if($cut[59] =~ /\((\d+)\/(\d+)\)\/(\d+)\)/){
			$pep3_num1=$1;
                        $pep3_num2=$3;
			$pep3_num3=$2;
			$pep3_rate=sprintf("%.2f", $pep3_num1/$pep3_num2);
			if($pep3_rate>=0.6  &&  $pep3_num2 <=7 && $pep3_num2>=5){
				$pep3_flag=1;
			}elsif($pep3_rate>=0.4  &&  $pep3_num2>=8 && $pep3_num1>$pep3_num3){
				$pep3_flag=1;
			}
		}elsif($cut[59] =~ /\((\d+)\/(\d+)\)/ )
		{
			$pep3_num1=$1;
			$pep3_num2=$2;
			$pep3_rate=sprintf("%.2f", $pep3_num1/$pep3_num2);
			if($pep3_rate>=0.6  &&  $pep3_num2 <=7 && $pep3_num2>=5){
                                $pep3_flag=1;
                        }elsif($pep3_rate>=0.4  &&  $pep3_num2>=8 ){
                                $pep3_flag=1;
                        }
			#$pep3_flag=1 if($pep3_num2>=7 && $pep3_rate>=0.6);
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
			elsif($pep3=~/c\.(\d+)del/  || $pep3=~/c\.(\d+)_(\d+)del/ || $pep3=~/c\.(\d+)_(\d+)ins/  || $pep3=~/c\.(\d+)dup/ ) #sunxq 190521
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
		#		$pep3_flag++ if($pep3_num1>0);
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
		if( $pep3_flag>0 ){ 
			$flag_PVS = 1; 
		}else{
			push @acmg,"PM4";
		}
		#if($cut[8]=~/^frameshift/ && $pep3_flag==0 ){push @acmg,"PM4";}
	}

	if($flag_PVS==1){ push @acmg,"PVS"; }
	
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
	
	#my $cont_ditacion_term=$control_ditaction_ratio{ join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-',1256); #NOVO1256样本检出率
	#my $case_ditacion_term=$case_ditacion_ratio{     join"\t",@cut[0,1,3,4]} // join"\t",(0,0,'-','-');
	
	if( ( ( (split/\t/,$case_ditacion_term)[1]==1 && (split/\t/,$case_ditacion_term)[2] eq $sample ) || (split/\t/,$case_ditacion_term)[1]==0 )  && (split/\t/,$cont_ditacion_term)[1]==0 )
	{
		$flag_PM2=1 if( $cut[44] eq "snp");
	}
	
	if($flag_PM2==1){ push @acmg,"PM2" if( (split /\|/, $cut[42])[1] eq "-" ); }
	
	my $ratio1=$cut[19];
	$ratio1=0 if($cut[19] eq "." && $db_dep_flag1 eq "Y" && $db_dep_12>=20 );
	my $ratio2=$cut[21];
	$ratio2=0 if($cut[21] eq "." && $db_dep_flag2 eq "Y" && $db_dep_22>=20 );
	my $ratio3=$cut[20];
	$ratio3=0 if($cut[20] eq ".");
	my $ratio4=$cut[22];
	$ratio4=0 if($cut[22] eq ".");
	
	my $ratio_max = max($ratio1,$ratio2);	
	#if(!exists $disease_ratio{$disease}){print "$disease\t$site_info\n";}
	#if($ratio_max < sqrt($disease_ratio{$disease}) && $adar=~/AR/i && $adar!~/AD/i )
	#{
	#	push @acmg,"PM2";
	#}

	my $ratio_min;
	if( $ratio1 ne "." && $ratio2 ne "." )
	{
		$ratio_min = min($ratio1,$ratio2,$ratio3,$ratio4);
		push @acmg,"BS1" if($ratio_min>=0.001); 
		
		#if($ratio_min > ( 1-sqrt(1 - $disease_ratio{$disease}) ) && $adar=~/AD/i )
		#{
		#	push @acmg,"BS1";
		#}

		#if($ratio_min > sqrt($disease_ratio{$disease}) && $adar=~/AR/i && $adar!~/AD/i )
		#{
		#	push @acmg,"BS1";
		#}
	}
	return @acmg;
}
#PVS，PS，PM，PP，BS ，BP，BA，
sub PathogenicityJudge{
	my($evidatabase,$eviacmg)=@_;
	my $tagDab=$1 if($evidatabase=~m/-(\S)$/);
	#print "$evidatabase\t$tagDab\t$eviacmg\n";
	$evidatabase=~s/-(\S)$//;
	my %evidenceCluase=();
	my $PathogenicityResult="无";
	my $flag=0;
	my $evicluasepren='';
	my $pathnum=0;
	my $benign=0;
	if($evidatabase !~/null/i ){# && ($tagDab eq 'N' || $tagDab eq 'R') ){
			my @subevi=split /\+/,$evidatabase;
			foreach my $vidvien (@subevi) {
				$evicluasepren.="$vidvien+";
				if($vidvien =~/致病/){
					$PathogenicityResult=$vidvien;
					$flag=1;
				}else{
					$vidvien=~s/\d$//;
					$pathnum++ if($vidvien=~/^P/);
					$benign++ if($vidvien=~/^B/);
					$evidenceCluase{$vidvien}++;
				}
			}
	}
	my $n=@$eviacmg;
	if($n>0){
		foreach my $subevi (@$eviacmg) {
				$evicluasepren.="$subevi+";
				$subevi=~s/\d$//;
				$pathnum++ if($subevi=~/^P/);
				$benign++ if($subevi=~/^B/);
				$evidenceCluase{$subevi}++;
		}
		#$evicluasepren=~s/\+$//;
	}
	$evicluasepren=~s/\+$//;
	if($flag==0){
		if($pathnum>0 && $benign==0){
			my ($pvsnum,$psnum,$pmnum,$ppnum)=(0,0,0,0);
			if(exists $evidenceCluase{PVS}){$pvsnum=$evidenceCluase{PVS};}
			if(exists $evidenceCluase{PS} ){$psnum=$evidenceCluase{PS};}
			if(exists $evidenceCluase{PM} ){$pmnum=$evidenceCluase{PM};}
			if(exists $evidenceCluase{PP} ){$ppnum=$evidenceCluase{PP};}
			if($pvsnum>0){
				if($psnum>0){
					$PathogenicityResult="致病" ;
				}else{
					if($pmnum>=2){
						$PathogenicityResult="致病" ;
					}elsif($pmnum==1){
						if($ppnum>0){
							$PathogenicityResult="致病" ;
						}else{
							$PathogenicityResult="可能致病" ;
						}
					}else{
						if($ppnum>=2){
							$PathogenicityResult="致病" ;
						}else{
							$PathogenicityResult="临床意义未明1级" ;
						}
					}
				}
			}
			if($pvsnum==0 && $psnum>0){
				if($psnum>=2){
					$PathogenicityResult="致病" ;
				}else{
					if($pmnum>=3){
						$PathogenicityResult="致病" ;
					}elsif($pmnum==2){
						if($ppnum>=2){
							$PathogenicityResult="致病" ;
						}else{
							$PathogenicityResult="可能致病" ;
						}
					}elsif($pmnum==1){
						if($ppnum>=4){
							$PathogenicityResult="致病" ;
						}else{
							$PathogenicityResult="可能致病" ;
						}
					}else{
						if($ppnum>=2){
							$PathogenicityResult="可能致病" ;
						}else{
							$PathogenicityResult="临床意义未明1级" ;
						}
					}
				}
			}
			if($pvsnum==0 && $psnum==0 && $pmnum>0){
				if( $pmnum>=3){
					$PathogenicityResult="可能致病" ;
				}elsif($pmnum==2){
					if($ppnum>=2){
						$PathogenicityResult="可能致病" ;
					}else{
						$PathogenicityResult="临床意义未明1级" ;
					}
				}elsif($pmnum==1){
					if($ppnum>=4){
						$PathogenicityResult="可能致病" ;
					}elsif($ppnum>=2){
						$PathogenicityResult="临床意义未明1级" ;
					}else{
						$PathogenicityResult="临床意义未明2级" ;
					}
				}
			}
			if($pvsnum==0 && $psnum==0 && $pmnum==0 && $ppnum>0){
				if($ppnum>=3){
					$PathogenicityResult="临床意义未明2级" ;
				}else{
					$PathogenicityResult="临床意义未明3级" ;
				}
			}
		}elsif($pathnum==0 && $benign>0){
			my ($banum,$bsnum,$bpnum)=(0,0,0);
			if(exists $evidenceCluase{BA} ){$banum=$evidenceCluase{BA};}
			if(exists $evidenceCluase{BS} ){$bsnum=$evidenceCluase{BS};}
			if(exists $evidenceCluase{BP} ){$bpnum=$evidenceCluase{BP};}
			if($banum>0){
				$PathogenicityResult="良性" ;
			}
			if($banum==0 && $bsnum>0 ){
				if($bsnum >=2){
					$PathogenicityResult="良性" ;
				}elsif($bsnum ==1){
					if($bpnum>=1){
						$PathogenicityResult="可能良性" ;
					}else{
						$PathogenicityResult="临床意义未明5级" ;
					}
				}
			}
			if($banum==0 && $bsnum==0 && $bpnum >0 ){
				if( $bpnum >=2){
					$PathogenicityResult="可能良性" ;
				}else{
					$PathogenicityResult="临床意义未明5级" ;
				}
			}
		}elsif($pathnum>0 && $benign>0){
			$PathogenicityResult="临床意义未明4级" ;
		}
	}
	if($evicluasepren eq ""){
		$evicluasepren="无";
		$PathogenicityResult="临床意义未明4级" ;
	}
	return $PathogenicityResult,$evicluasepren;
}
