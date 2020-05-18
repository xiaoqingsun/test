use strict;
use warnings;
use utf8;
use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Encode;
use Spreadsheet::XLSX;
use List::Util qw/max min sum maxstr minstr shuffle/;

my ($out,$types,$cfg);
GetOptions(
        'out=s' => \$out,
        'type=s' => \$types,
        'cfg=s' => \$cfg,
);


my $usage=<<END;

perl $0 
        -path   需要检测样本的路径list文件，例如/ifs/TJPROJ3/HEALTH/Project/2C/160926_ST-E00192_0305_BH35GVALXX/cr/NHE0000000649/NHE0000000649/04YC/NHE0000000649.final
        -type   产品类型,YC | YC133(遗传病)、GXY、CS、CR 、NNY
        -out    输出路径,输出相应产品的true的点，且YC，GXY，CS会与最新位点数据库匹配,返回与数据库匹配的结果
END

my $typee=$types;
$types="IDT";
$typee=~s/IDT/WES/g;

my @type2=();
my %sampledata=();
push @type2,"GXYRA" if($typee =~/GXY/);
push @type2,"CSRA" if($typee =~/CS/); 
push @type2,"DMRA" if($typee =~/DM/);
push @type2,"WESRA" if($typee =~/WES/);

foreach my $type (@type2){
	my $base=judge_lastest_sitedb($type);
	if($type eq "GXYRA"){
		read_samdatabase_gxy($base);
	}else{
		read_samdatabase($base);
	}
}

##################################
open RT,"<:utf8", "$out/$types.right.txt" || die "can't open the $out/$types.right.txt!";
open PD, "> $out/$types.pdf.txt" || die "can't open the $out/$types.pdf.txt!";
my $head=<RT>; chomp $head; 
my @h=split /\t/,$head; my $hn=@h; $hn=$hn-1;
my $headnew=join("\t",@h[1..$hn-2]);
print PD "#$headnew\t解读人\ttime\n";
my %sample_flag;
my %sample_id;
while(<RT>){
	chomp;
	my @t=split /\t/;
	next if(/^$/);
	next if(/^#/);
        $sample_id{s}{$t[1]}=1;
        if($t[5] eq '.'){$t[5]='-';}
        my $key=join("_",@t[1..3,5,6,84]);
        $sample_id{k}{$key}=1;
	if(exists $sampledata{$key}{key}){ 
		$sample_flag{$t[1]}{k}=1; 
		if($sampledata{$key}{key} =~/report/i || $sampledata{$key}{key} =~/zhubiao/i ){ #||  $sampledata{$key}{key} =~/negative/i ){
			print PD "$sampledata{$key}{info}\n";
			$sample_flag{$t[1]}{f}=1;  
		}elsif($sampledata{$key}{key} =~/del/i){
			print PD "####$sampledata{$key}{info}\n";
		}elsif($sampledata{$key}{key} =~/beizhu/i){
			print PD "###$sampledata{$key}{info}\n";
		}elsif($sampledata{$key}{key} =~/fubiao/i){
			print PD "##$sampledata{$key}{info}\n";
		}elsif($sampledata{$key}{key} =~/negative/i){
			next;
		}else{
			die "请核实$key\备注为$sampledata{$key}{key}\n";
		}
	}else{
		shift @t;
		#print $t[-1];
		my $linenew=join("\t",@t);
		print PD "N$linenew\n";
	}
}
close RT;
foreach my $key(keys %sampledata){
	my @a=split /_/,$key;
        #print "##$key\n";
	if(exists $sample_id{s}{$a[0]} && !exists $sample_id{k}{$key}){
		$sample_flag{$a[0]}{k}=1;
		if($sampledata{$key}{key} =~/report/i || $sampledata{$key}{key} =~/zhubiao/i ){# || $sampledata{$key}{key} =~/negative/i ){
                        print PD "$sampledata{$key}{info}\n";
                        $sample_flag{$a[0]}{f}=1;  print "###$a[0]\n";
                }elsif($sampledata{$key}{key} =~/del/i){
                        print PD "####$sampledata{$key}{info}\n";
                }elsif($sampledata{$key}{key} =~/beizhu/i){
                        print PD "###$sampledata{$key}{info}\n";
                }elsif($sampledata{$key}{key} =~/fubiao/i){
                        print PD "##$sampledata{$key}{info}\n";
		}elsif($sampledata{$key}{key} =~/negative/i){
			next;
                }else{
                        die "请核实$key\备注为$sampledata{$key}{key}\n";
                }	
	}
}
foreach my $key(keys %sample_flag){
	my $key2=$key."_-_-_-_-_-";
	if(!exists $sample_flag{$key}{f} && exists $sample_flag{$key}{k}){
		if(exists $sampledata{$key2}{key} && $sampledata{$key2}{key} =~/negative/i){
			print $key;
			print PD "$sampledata{$key2}{info}\n";
		}else{
		my $null="-\t" x 92;
		print PD "$key\t$null\n";
		}
	}
}
close PD;

####################################
sub read_samdatabase{
	my $shuju=shift;
	#my %hash;
        open IN, "<:utf8", "$shuju" || die "can't open the file";
	while(<IN>){
                chomp;
                my @tmp=split/\t/,$_;
		my $n=@tmp;
		$n=$n-1;
		if($tmp[6] eq '.'){$tmp[6]='-';}
                my $keysam=join("_",@tmp[2..4,6,7,85]);
		if(!exists $sampledata{$keysam}){
        		$sampledata{$keysam}{info}=join("\t",@tmp[2..$n]);
			$sampledata{$keysam}{key}=$tmp[0];  
		}
	}
	close IN;
	#return %hash;
}

sub  read_samdatabase_gxy{
	my $shuju=shift;
	open IN, "<:utf8", "$shuju" || die "can't open the file";
	while(<IN>){
		chomp;
		my @tmp=split/\t/,$_;
		my $n=@tmp;
		if($tmp[6] eq '.'){$tmp[6]='-';}
		if($tmp[7] eq '.'){$tmp[7]='-';}
		my $keysam=join("_",@tmp[2..4,6,7,84]);
		my @site_info=(@tmp[2..71],$tmp[77],@tmp[72..77],@tmp[78..88],$tmp[91],@tmp[105..107],$tmp[$n-2],$tmp[$n-1]);
		if(!exists $sampledata{$keysam}){
			$sampledata{$keysam}{info}=join("\t",@site_info);
			$sampledata{$keysam}{key}=$tmp[0];
		}
	}
	close IN;
}

sub judge_lastest_sitedb{
        my $product_type=shift;
        my %hash_product_sitepath=();
        my @all_date;
        my $lastest_date;
        my $lastest_sitedb_path;
	my $gxy_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/GXY_report_base/";
	my $cusi_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/CuS_report_base/";
        my $dm_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/DM_report_base/";
	my $wes_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_report_base/";

	$hash_product_sitepath{'GXYRA'}=$gxy_report_all_path;
        $hash_product_sitepath{'CSRA'}=$cusi_report_all_path;
        $hash_product_sitepath{'DMRA'}=$dm_report_all_path;
	$hash_product_sitepath{'WESRA'}=$wes_report_all_path;
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


