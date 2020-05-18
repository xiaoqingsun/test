#!/usr/bin/perl -w
use strict;
use utf8;
use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use File::Path qw(make_path);
use List::Util qw/max min sum maxstr minstr shuffle/;


################ database path##
my $gxy_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY";
my $gxy_report_sample_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/GXY_report/";
my $gxy_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/GXY_report_base/";
my $gxy_trace="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY_trace/";
my $gxy_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_GXY_fubiao";
my $dm_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/DM/";
my $dm_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/DM_fubiao/";
my $dm_report_sample_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/DM_report/";
my $dm_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/DM_report_base/";
my $dm_trace="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_DM_trace/";
my $cusi_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_CuS";
my $cusi_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_CuS_fubiao";
my $cusi_report_sample_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/CuS_report/";
my $cusi_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/CuS_report_base/";
my $cusi_trace="/PUBLIC/pipline/database/knowledge_database/sitedatabase/Xinkang_CuS_trace/";
my $wes_sitedb_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES/";
my $wes_sitedb_path_fubiao="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_fubiao/";
my $wes_report_sample_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_report/";
my $wes_report_all_path="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_report_base/";
my $wes_trace="/PUBLIC/pipline/database/knowledge_database/sitedatabase/WES_trace/";
my $genetic_black_box="/NJPROJ1/HEALTH/database/knowledge_database/sitedatabase/genetic_heihe";
my $csgene_disease_file="/PUBLIC/pipline/database/siteFilterDB/product/CuS/XK_cus_disease_subtype.info";
my $gxygene_disease_file="/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/gene_deseases/170718_v4_gene_desease";
my $dmgene_disease_file="/PUBLIC/work/sunxiaoqing/cs_report/DM_disease_list";
#########################

########################read database##
my @producttype=("GXY","CS","DM","WES");
foreach my $type (@producttype){
	my (%database,%sample_data,%all_disease);
	my $falg1=$type."RS"; my $falg2=$type."RA"; my $falg3=$type."T"; my $falg4=$type."FB";
	my $run_flag=1;

	my ($sample_report_database,$uptime_sa)=judge_lastest_sitedb($falg1);
	my ($all_report_database,$uptime_saAll)=judge_lastest_sitedb($falg2);
	my ($readdatabase,$uptime)=judge_lastest_sitedb("$type");
	my ($readtrace,$uptimetrace)=judge_lastest_sitedb($falg3);
	my ($readfubao,$uptime_fubiao)=judge_lastest_sitedb($falg4);
	if($uptime_sa   &&  $uptime  && $uptime_sa <=$uptime){
		print  "位点解读库日期晚于更新样本库日期，$type 库将不会更新！\n"; $run_flag=0; 
	}
	if($uptime_sa   &&  $uptimetrace  && $uptime_sa <= $uptimetrace){
		print  "追踪库日期晚于更新样本位点库日期，$type 库将不会更新！\n"; $run_flag=0; 
	}
	if($uptime_sa   &&  $uptime_saAll  && $uptime_sa <=$uptime_saAll){
		print "所有样本出报告位点库日期晚于更新样本库日期，$type 库将不会更新！\n"; $run_flag=0; 
	}

#	print "$sample_report_database,$uptime_sa\n$all_report_database,$uptime_saAll\n$readdatabase,$uptime\n$readtrace,$uptimetrace\n";

        next if($run_flag==0);

	if($type eq 'GXY'){
        	%all_disease=&gxyread_disease($gxygene_disease_file); 
		%sample_data=readGXY_sample_report_data($sample_report_database,\%all_disease,$uptime_sa);
        	%database=&updata_read_base($readdatabase,$readfubao,$gxy_sitedb_path,84,'XK_gxybase',\%sample_data,$uptime_sa);
	        &updata_read_base_fubiao(\%database,$gxy_sitedb_path_fubiao,84,'XK_gxybase_fubiao',\%sample_data,$uptime_sa,$readfubao);
        	&updata_readtrace_base($readtrace,$gxy_trace,'gxytracebase',\%sample_data,$uptime_sa);
	        &updata_all_report_database($all_report_database,$gxy_report_all_path,'gxy_reportAllbase',\%sample_data,$uptime_sa);
	}elsif($type eq 'CS'){
        	%all_disease=&read_disease($csgene_disease_file);
		%sample_data=readCS_sample_report_data($sample_report_database,\%all_disease,$uptime_sa,$type);
		%database=&updata_read_base($readdatabase,$readfubao,$cusi_sitedb_path,85,'XK_cusbase',\%sample_data,$uptime_sa);
        	&updata_read_base_fubiao(\%database,$cusi_sitedb_path_fubiao,85,'XK_CSbase_fubiao',\%sample_data,$uptime_sa,$readfubao);
		&updata_readtrace_base($readtrace,$cusi_trace,'custracebase',\%sample_data,$uptime_sa);
		&updata_all_report_database($all_report_database,$cusi_report_all_path,'CuS_reportAllbase',\%sample_data,$uptime_sa);
	}elsif($type eq 'DM'){
        	%all_disease=&gxyread_disease($dmgene_disease_file);
	        %sample_data=&readCS_sample_report_data($sample_report_database,\%all_disease,$uptime_sa,$type);
        	%database=&updata_read_base($readdatabase,$readfubao,$dm_sitedb_path,85,'XK_dmbase',\%sample_data,$uptime_sa);
	        &updata_read_base_fubiao(\%database,$dm_sitedb_path_fubiao,85,'XK_DMbase_fubiao',\%sample_data,$uptime_sa,$readfubao);
       		&updata_readtrace_base($readtrace,$dm_trace,'dmtracebase',\%sample_data,$uptime_sa);
	        &updata_all_report_database($all_report_database,$dm_report_all_path,'dm_reportAllbase',\%sample_data,$uptime_sa);
	}elsif($type eq 'WES'){
		%sample_data=&readCS_sample_report_data($sample_report_database,\%all_disease,$uptime_sa,$type);
		%database=&updata_read_base($readdatabase,$readfubao,$wes_sitedb_path,85,'XK_wesbase',\%sample_data,$uptime_sa);
		&updata_read_base_fubiao(\%database,$wes_sitedb_path_fubiao,85,'XK_WESbase_fubiao',\%sample_data,$uptime_sa,$readfubao);
		&updata_readtrace_base($readtrace,$wes_trace,'westracebase',\%sample_data,$uptime_sa);
		&updata_all_report_database($all_report_database,$wes_report_all_path,'wes_reportAllbase',\%sample_data,$uptime_sa);
	}
}
##################################
sub gxyread_disease{
	my ($file)=@_;
        my %hash;
        open(IN,"<:utf8",$file)||die;
        while(<IN>){
                chomp;
                my @t=split /\t/;
                $hash{"$t[0]\t$t[1]"}=1;
        }
        close IN;
        return %hash;
}
sub read_disease{
	my ($file)=@_;
	my %hash;
	open(IN,"<:utf8",$file)||die;
	while(<IN>){
		chomp;
		my @t=split /\t/;
		$hash{"$t[1]\t$t[3]"}=1;
	}
	close IN;
	return %hash;
}
sub updata_all_report_database{
	my ($basefile,$baspath,$basn,$sample_datas,$uptime_saa)=@_;
        make_path("$baspath/$uptime_saa");
	my $basefilename="$baspath/$uptime_saa/$basn$uptime_saa.right.txt";
	open(OUT,">$basefilename")||die;
	my %hash3;
	if($basefile){
		open(IN,"<:utf8",$basefile)||die;
		while(<IN>){
                        chomp;
			my @tmp=split /\t/,$_;
			if(!exists $$sample_datas{samsite}{$tmp[2]}){
				$hash3{$tmp[2]}=1;
                        	print OUT "$_\n"
			}
		}
		close IN;
	}
	foreach my $key(keys %{$$sample_datas{samsite}}){
		print OUT "$$sample_datas{samsite}{$key}"; 
	}
	close OUT;
}

sub updata_readtrace_base{
	my ($basefile,$baspath,$basna,$sample_datas,$uptime_saa)=@_;
	make_path("$baspath/$uptime_saa"); 
	my $basefilename="$baspath/$uptime_saa/$basna$uptime_saa.right.txt";	
	open(OUT,">$basefilename")||die;
	if($basefile){ 
		open(IN,"<:utf8",$basefile)||die;
		while(<IN>){
			chomp;
			print OUT "$_\n";
		}
		close IN;
	}
	foreach my $key(keys %{$$sample_datas{reportsite}}){
		print OUT "$$sample_datas{reportsite}{$key}\n";
	}
	close OUT;
}

sub updata_read_base_fubiao{
	my ($hash,$baspath,$discol,$basn,$sample_datas,$uptime_saa,$readfubaoo)=@_;
	make_path("$baspath/$uptime_saa");
        my $basefilename="$baspath/$uptime_saa/$basn$uptime_saa.right.txt";
        open(OUT,">$basefilename")||die;
        my %hash4;
        if($readfubaoo){
        	open(IN,"<:utf8",$readfubaoo)||die;
		while(<IN>){
			chomp;
			my @tmp=split /\t/,$_;
                	my $key1=join("_",@tmp[3,4,6,7,$discol]);
			my $key2=join("_",@tmp[3,4,6,7]);
                	if(!exists $$hash{$key2} || $_=~/Flag/){  ### panqi 主表无此位点
                        	$hash4{$key2}=1;
				if(exists $$sample_datas{fubiaosite}{$key1}){        ### panqi 190106
                                	print OUT "T1$$sample_datas{fubiaosite}{$key1}\n";
                        	}else{
					print OUT  "$_\n";
				}
			}
		}
		close IN;
	}
	foreach my $key(keys %{$$sample_datas{fubiao}}){
		if(!exists $$hash{$key}  && !exists $hash4{$key}){  ### panqi 主表附表都无此位点
			print OUT $$sample_datas{fubiao}{$key};
		}
	}
	close OUT;
}

sub updata_read_base{
	my ($basefile,$readfubaoo,$baspath,$discol,$basn,$sample_datas,$uptime_saa)=@_;
	make_path("$baspath/$uptime_saa");
	my $basefilename="$baspath/$uptime_saa/$basn$uptime_saa.right.txt";
	my %hash2;
	my %hash3;
	my %hash4;
	open(OUT,">$basefilename")||die;  
	if($readfubaoo){
		open(IN,"<:utf8",$readfubaoo)||die;
		while(<IN>){
			chomp;
			my @tmp=split /\t/,$_;
			my $key=join("_",@tmp[3,4,6,7,$discol]);
			my $key2=join("_",@tmp[3,4,6,7]);
			$hash4{$key2}.=join("\t",@tmp[0..89])."\n" if($discol == 85 && !exists $$sample_datas{site}{$key} );
			$hash4{$key2}.=join("\t",@tmp[0..119])."\n" if($discol == 84 && !exists $$sample_datas{site}{$key} );
		}
		close IN;
	}
	if($basefile){
		open(IN,"<:utf8",$basefile)||die;
		while(<IN>){
			chomp;
			my @tmp=split /\t/,$_;
			my $key=join("_",@tmp[3,4,6,7,$discol]);
			my $key2=join("_",@tmp[3,4,6,7]);
			$hash2{$key}=1;
			$hash3{$key2}=1;
			if(exists $$sample_datas{site}{$key}){
				print OUT "T1$$sample_datas{site}{$key}\n";
			}else{
				print OUT "$_\n";
			}
		}
		close IN;
	}
	
	foreach my $key(keys %{$$sample_datas{site}}){
		my @sss=split /_/,$key;
		my $keysub=(join "_",@sss[0,1,3,4]);
		next if(exists $hash3{$keysub});
		$hash3{$keysub}=1;
		if (!exists $hash2{$key})
		{
			print OUT "Z$$sample_datas{site}{$key}\n";
			print OUT $hash4{$keysub} if(exists $hash4{$keysub});
		}
	}
	close OUT;
	return %hash3;
}

sub readGXY_sample_report_data{
	my ($shuju,$all_diseases,$uptime_saa)=@_;#shift; print "#$shuju\n";
        open IN, "<:utf8", "$shuju" || die "can't open the file";
        my %hash;
        my %head;
        my $aa=<IN>;chomp $aa;
        my @heads=split /\t/,$aa;
        for(my $i=1;$i<@heads;$i++){
                $head{$heads[$i]}=$i;
        }
        chomp $aa;
	while(<IN>){
		chomp;
                my @tmp=split/\t/,$_;
                my $gene_disea="$tmp[86]\t$tmp[84]";  
                die "请检查$gene_disea 是否和数据库一致 " if($tmp[1] !~/Flag/i  && !exists $$all_diseases{$gene_disea} &&  $tmp[0] !~/Negative/i);
                my ($key,$keysam);
                $key=join("_",@tmp[3,4,6,7,84]);
                $keysam=join("_",@tmp[2..4,6,7,84]);
                $hash{samsite}{$tmp[2]}.="$_\t$uptime_saa\n";  
                if($tmp[0] =~/report/i){ 
                        my $person=$head{'解读人'};
                        $hash{site}{$key}=join("\t",$tmp[$person],@tmp[1..119]);
                        $tmp[0]=$uptime_saa;
                        $hash{reportsite}{$keysam}=join("\t",@tmp);
                }elsif($tmp[0] =~/fubiao/i){ 
                        my $key4=join("_",@tmp[3,4,6,7]);
			$hash{fubiao}{$key4}.="$_\t$uptime_saa\n";
			$hash{fubiaosite}{$key}="$_\t$uptime_saa\n";  ### panqi 190106
		}
	}
	close IN;
	return %hash;
}
sub readCS_sample_report_data{
	my ($shuju,$all_diseases,$uptime_saa,$typee)=@_;#my $shuju=shift;
	open IN, "<:utf8", "$shuju" || die "can't open the file";
	my %hash;
	my %head;
	my $aa=<IN>;chomp $aa;
	my @heads=split /\t/,$aa;
	for(my $i=1;$i<@heads;$i++){
		$head{$heads[$i]}=$i;
	}
        chomp $aa;
	while(<IN>){
                chomp;
                my @tmp=split/\t/,$_;
		my $gene_disea="$tmp[87]\t$tmp[85]";
		die "请检查$gene_disea 是否和数据库一致 " if($typee ne 'WES'&& $tmp[1] !~/product/i  && !exists $$all_diseases{$gene_disea} && $tmp[0] !~/Negative/i);
                my ($key,$keysam);
		$key=join("_",@tmp[3,4,6,7,85]);
                $keysam=join("_",@tmp[2..4,6,7,85]);
		$hash{samsite}{$tmp[2]}.="$_\t$uptime_saa\n";
		if($tmp[0] =~/report/i){
			my $person=$head{'解读人'};
			$hash{site}{$key}=join("\t",$tmp[$person],@tmp[1..89]); 
			$tmp[0]=$uptime_saa;
			$hash{reportsite}{$keysam}=join("\t",@tmp);
		}elsif($tmp[0] =~/fubiao/i){
			my $key4=join("_",@tmp[3,4,6,7]);
			$hash{fubiao}{$key4}.="$_\t$uptime_saa\n";
			$hash{fubiaosite}{$key}="$_\t$uptime_saa\n"; ### panqi 190106                
		}
	}
	close IN;
        return %hash;
}


sub judge_lastest_sitedb{
        my $product_type=shift;
        my %hash_product_sitepath=();
        my @all_date;
        my $lastest_date;
        my $lastest_sitedb_path;
        $hash_product_sitepath{'GXY'}=$gxy_sitedb_path;
        $hash_product_sitepath{'GXYRS'}=$gxy_report_sample_path;
        $hash_product_sitepath{'GXYRA'}=$gxy_report_all_path;
        $hash_product_sitepath{'GXYT'}=$gxy_trace;
        $hash_product_sitepath{'DM'}=$dm_sitedb_path;
        $hash_product_sitepath{'DMRS'}=$dm_report_sample_path;
        $hash_product_sitepath{'DMRA'}=$dm_report_all_path;
        $hash_product_sitepath{'DMT'}=$dm_trace;
        $hash_product_sitepath{'CS'}=$cusi_sitedb_path;
        $hash_product_sitepath{'CSRS'}=$cusi_report_sample_path;
        $hash_product_sitepath{'CSRA'}=$cusi_report_all_path;
	$hash_product_sitepath{'CST'}=$cusi_trace;
        $hash_product_sitepath{'WES'}=$wes_sitedb_path;
        $hash_product_sitepath{'WESRS'}=$wes_report_sample_path;
        $hash_product_sitepath{'WESRA'}=$wes_report_all_path;
        $hash_product_sitepath{'WEST'}=$wes_trace;
        $hash_product_sitepath{'black'}=$genetic_black_box;
        #$hash_product_sitepath{'fubiao'}=$cusi_sitedb_path_fubiao;
        $hash_product_sitepath{'CSFB'}=$cusi_sitedb_path_fubiao;
        $hash_product_sitepath{'GXYFB'}=$gxy_sitedb_path_fubiao;
        $hash_product_sitepath{'DMFB'}=$dm_sitedb_path_fubiao;
	$hash_product_sitepath{'WESFB'}=$wes_sitedb_path_fubiao;
	if( exists $hash_product_sitepath{$product_type}){
                my @all_path=glob "$hash_product_sitepath{$product_type}/*";
                foreach(@all_path){
                        my$date=(split/\//)[-1];
                        if($date=~/^\d+$/i){
                                push @all_date,$date;
                        }
                }
                $lastest_date=max(@all_date);
		#if($product_type eq 'CS'){ $uptime=$product_type; }
                my @lastest_sitedb_path=glob "$hash_product_sitepath{$product_type}/$lastest_date/*right.txt";
                foreach(@lastest_sitedb_path){
                        if(/right/i){
                                $lastest_sitedb_path=$_;
                        }
                }
        }
        print "$lastest_sitedb_path\n";
        print "$product_type\n";
        return $lastest_sitedb_path,$lastest_date;
}

