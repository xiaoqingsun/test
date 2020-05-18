#!/usr/bin/perl -w
use strict;

use warnings;
use Encode;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use feature qw/say/;
sub printerr;

my (%psudePos,%HSPos,%pLI_oe_gene,%repeat_site,%variant,%variant2,%gnomdepth,%exacdepth);
my ($input,$bedfile,$type,$bed_file,$posarry,@posarrys);

GetOptions(
    'in=s' => \$input,
    'bed=s' => \$bedfile,
    'pos=s' => \$posarry
#    'outdir=s' => \$outdir,
);

die `pod2text $0` if ($input eq '');
my $file_psude ||= "/PUBLIC/pipline/database/clingen_Haploinsufficiency/file.psudegene.pos";
my $file_pLI_oe ||= "/PUBLIC/pipline/database/clingen_Haploinsufficiency/file.pLI_oe.txt";
my $file_HS ||= "/PUBLIC/pipline/database/clingen_Haploinsufficiency/file.HS";
my $dbEXAC ||= "/PUBLIC/work/sunxiaoqing/EXAC/Panel.all.coverage.txt";
my $dbexomes ||= "/PUBLIC/work/sunxiaoqing/gnomad/gnomad.exomes.r2.0.1.coverage.txt";
die "数据库没有idx文件 " unless (-f "$dbEXAC.idx" && -f "$dbexomes.idx");

if($bedfile=~/IDT/){
	$bed_file="/PUBLIC/pipline/database/bed/idt.repeatbase.bed";
}else{
	$bed_file="/PUBLIC/pipline/database/bed/panel.repeatbase.bed";
}


if($input=~/INDEL/){
        $type="indel";
}elsif($input=~/SNP/i){
	$type="snp";
}else{
	$type="";
}


open(IN,$input) or die;
while(<IN>){
	chomp;
	next if(/^Chr/);
	my @t=split /\t/;
	my ($chr, $start, $end) = @t[0..2];
	$variant2{$chr, $start} = 1;
	$variant2{$chr, $end} = 1;
	for($start-10..$end+10){
		$variant{$chr, $_} = 1;
	}
}
close IN;


`mv $input $input.bak`;
if($input=~/INDEL/){
	$type="indel";
	&read_repeatsite;
}
&read_psude_region;
&read_HS_region;
&read_pLI_oe_gene;
%gnomdepth=&read_depth_db($dbexomes,\%variant,\%variant2,$type);
%exacdepth=&read_depth_db($dbEXAC,\%variant,\%variant2,$type);

@posarrys=split /,/,$posarry;
open IN, "$input.bak" or die $!;
open OUT,">","$input" or die $!;
while(<IN>){
	chomp;	
	my @t=split /\t/;
	my ($chr, $start, $end) = @t[0..2];
	if(/^CHR/){
                
		$t[$posarrys[0]]="Pseudgene";
		$t[$posarrys[1]]="HS_pLI_oe";
		$t[$posarrys[2]]="RepeatBase";
		$t[$posarrys[3]]="DBdepth";
	}else{
		my $peudeinfo="";
		my $HSinfo="";
		my $pLI_oe_info="";
		my $repeat="";
		if(exists $psudePos{$t[0]}){
			for(my $i=0;$i<@{$psudePos{$t[0]}};$i++){
				my $st=${$psudePos{$t[0]}}[$i][0];
				my $en=${$psudePos{$t[0]}}[$i][1];
				if($st>=$start && $st<=$end || $start>=$st && $start<=$en){
					$peudeinfo.="${$psudePos{$t[0]}}[$i][2];";
				}
			}
		}
		
		if(exists $HSPos{$t[0]}){
                        for(my $i=0;$i<@{$HSPos{$t[0]}};$i++){
                                my $st=${$HSPos{$t[0]}}[$i][0];
                                my $en=${$HSPos{$t[0]}}[$i][1];
                                if($st>=$start && $st<=$end || $start>=$st && $start<=$en){
                                        $HSinfo.="HS:${$HSPos{$t[0]}}[$i][2],${$HSPos{$t[0]}}[$i][3];";
                                }
                        }
                }
		if($type eq "indel"){
			if(exists $repeat_site{$t[0]}){
				for(my $i=0;$i<@{$repeat_site{$t[0]}};$i++){
					my $st=${$repeat_site{$t[0]}}[$i][0];
					my $en=${$repeat_site{$t[0]}}[$i][1];
					if($st>=$start && $st<=$end || $start>=$st && $start<=$en){
						$repeat.="repeat:".join("_",@{${$repeat_site{$t[0]}}[$i]});
					}
				}
			}
		}
		
		if(exists $pLI_oe_gene{$t[6]}){
			$pLI_oe_info=$pLI_oe_gene{$t[6]};
		}	

		$peudeinfo="-" if( $peudeinfo eq "");
		$HSinfo="-;" if($HSinfo eq "");
		$pLI_oe_info="-" if($pLI_oe_info eq "");
		$repeat="-" if($repeat eq "");

		my ($dep_gnomAD,$flag_cov_gnomAD)=&get_site_depth($chr, $start, $end,\%gnomdepth);
		my ($dep_Exac,$flag_cov_exac)=&get_site_depth($chr, $start, $end,\%exacdepth);
		my $depthLine="$dep_gnomAD//$dep_Exac//$flag_cov_gnomAD;$flag_cov_exac";

		$t[$posarrys[0]]=$peudeinfo;
		$t[$posarrys[1]]="$HSinfo$pLI_oe_info";
		$t[$posarrys[2]]=$repeat;	
		$t[$posarrys[3]]=$depthLine;	
	}
	my $line=join("\t",@t); 
	print OUT  "$line\n";
}
close IN;

########################sub roution
sub read_psude_region{
#PRKAR1A	NM_212471	11	17	66526418	66527949	atten83-86%相似
	open IN, $file_psude or die $!;
	while(<IN>){
		chomp;
		next if(/#/);
		my @t=split /\t/;
		my $inf="$t[0]:$t[1]:exon$t[2]:$t[6]";
		push @{$psudePos{$t[3]}},[$t[4],$t[5],$inf];
	}
	close IN;
}

sub read_HS_region{
#	chr1    10535002        10690815        PEX14   30
	open IN, $file_HS or die $!;
	while(<IN>){
                chomp;
		my @t=split /\t/;
		$t[0]=~s/chr//;
		push @{$HSPos{$t[0]}},[$t[1],$t[2],$t[3],$t[4]];	
	}
	close IN;
}

sub read_pLI_oe_gene{
# gene    pLI     oe_lof  oe_lof_lower    oe_lof_upper    gene_id HGNCid  HGNCsymblo      RefGene
# MED13   1.0000e+00      0.0000e+00      0.0000e+00      3.0000e-02      ENSG00000108510 HGNC:22474      MED13   MED13
	open IN,$file_pLI_oe or die $!;
	while(<IN>){
                chomp;
                my @t=split /\t/;
		my $info="pLI:$t[1];oe:$t[2]\[$t[3],$t[4]\]";
		my $gname=$t[0];
		$gname=$t[8] if($t[8] ne '-');
		$pLI_oe_gene{$gname}=$info;
	}
        close IN; 
}

sub read_repeatsite{
	open IN,$bed_file or die $!;	
	while(<IN>){
		chomp;
		my @t=split /\t/;
		push  @{$repeat_site{$t[0]}},[$t[4],$t[5],$t[6]];
	}
	close IN;	
}


#####################
sub get_site_depth{
	my ($chr, $start, $end,$hashDepth)=@_;
	my ($fronts,$varss,$behinds)=('','','');
	my @front=('-','-','-','-','-','-','-','-','-','-');
	my @behind=('-','-','-','-','-','-','-','-','-','-');
	my $length=$end-$start+1;
	my $var='-'x$length;
	my @vars=split //,$var;
	my $flag=0;
	my $i=0;
	for($start-10..$start-1){
		if(exists $$hashDepth{$chr,$_}){$front[$i]= $$hashDepth{$chr,$_};$flag++;}
		$i++;
	}
	
	$i=0;
	for($start..$end){
		if(exists $$hashDepth{$chr,$_}){$vars[$i]= $$hashDepth{$chr,$_};$flag++;}
		$i++;
	}

	$i=0;
	for($end+1..$end+10){
		if(exists $$hashDepth{$chr,$_}){$behind[$i]= $$hashDepth{$chr,$_};$flag++;}
		$i++;
	}

	$fronts=join(",",@front);
	$varss=join(",",@vars);
	$behinds=join(",",@behind);
	my $retrunDP="$fronts|$varss|$behinds";
	my $flagC="N";
	$flagC="Y" if($flag>0);
	return $retrunDP,$flagC;
}

sub read_depth_db{
	my($dbfile,$variant,$variant2,$class)=@_;
	my $flag_idx_search=0;
	my $indexfilter_threshold = 0.9;
	my ($BIN, $DBSIZE) = (0, 0);
	my %dbdepht;
	my $bb = {};			#a subset of %index, which corresponds to the input variants
	my %index = ();			#holds all bin index, used by DBM only and not by IDX 
	open(IDX, "$dbfile.idx")||die;
	my $line = <IDX>;
	if (not $line =~ m/BIN\t(\d+)\t(\d+)/) {
		printerr "WARNING: Malformed database index file $dbfile.idx.\n";
	} elsif ($2 != -s $dbfile) {                    #file version is different, do not use index file in this case
		printerr "WARNING: Your index file $dbfile.idx is out of date and will not be used. ANNOVAR can still generate correct results without index file.\n";
	} else {
		($BIN, $DBSIZE) = ($1, $2);
		foreach my $k ( keys %$variant2 ) {
			my ($chrom, $pos) = split ($;, $k);
			my $bin = $pos - ($pos % $BIN);
			$bb->{"$chrom\t$bin"} = 0;      #flag this bin to be searched later
		}

		my ($count_total, $count_search) = qw/0 0/;
		while ( $line = <IDX> ) {
			$line =~ s/[\r\n]+$//;
			my ( $chrom, $bin, $offset0, $offset1 ) = split (/\t/, $line);
			$chrom =~ s/^chr//;             #delete the chr in snp135, etc 
			defined $offset1 or next;       #invalid input line in the index file

			if (defined $bb->{"$chrom\t$bin"}) {
				$bb->{"$chrom\t$bin"} or $count_search++;
				$bb->{"$chrom\t$bin"} = "$offset0,$offset1";
			}
			$count_total++;
		}
		if ($count_total and $count_search/$count_total < $indexfilter_threshold) {
			$flag_idx_search++;
			printerr "NOTICE: Database index loaded. Total number of bins is $count_total and the number of bins to be scanned is $count_search\n";
		}
		if (not $flag_idx_search) {
                $bb = {1, join (',', 0, -s "$dbfile")};
                %index = (1, join (',', 0, -s "$dbfile"));
        	}

	        open (DB, $dbfile) or die "Error: cannot read from input database file $dbfile: $!\n";
	        printerr "NOTICE: Scanning filter database $dbfile...";

		foreach my $b (sort keys %$bb) {
			$bb->{$b} or next;					#DB does not have this bin to search against (value is zero)
			my ($chunk_min, $chunk_max) = split (/,/, $bb->{$b});
			$chunk_min-=2000;
			$chunk_max+=2000;	
			defined $chunk_max or die "Error: b=$b index=$index{$b} bb=$bb->{$b}\n";
			seek(DB, $chunk_min, 0);				#place file pointer to the chunk_min
			my $chunk_here = $chunk_min;
			<DB>;
			while (<DB>) {
				my $line=$_;
				my $line_length = length($_);			#calculate line length of the DB
				$chunk_here += $line_length;
				s/[\r\n]+$//;
				m/\S/ or next;					#skip empty lines in the database file (sometimes this occurs)
				m/^#/ and next;					#skip the comment line
				chomp;
				my @record = split (/\t/, $_);
				my ($chr, $start, $obs) = ($record[0],$record[1],$record[2]);
				$chr =~ s/^chr//i;		
				if($variant->{$chr, $start}){$dbdepht{$chr, $start}=$obs; }
				if ( $chunk_here > $chunk_max ) {
					last;
				}
			}
		}
	}
	return  %dbdepht;
}


sub printerr {
	print STDERR @_;
#	print LOG @_;
}

