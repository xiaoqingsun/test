use strict;
use warnings;
#use utf8;
#use Encode;
use autodie;
use Data::Dumper;
use File::Basename;
use Smart::Comments;
use Getopt::Long;
use Spreadsheet::XLSX;
use List::Util qw/max min sum maxstr minstr shuffle/;

my $cfg_file=shift;
my $version=shift;

my @cut;
my %id_hash;

open IN, $cfg_file or die $!;
while(<IN>)
{
	chomp;
	next if(/^样本编号/);
	@cut=split /\t/, $_;
	$id_hash{$cut[0]}=1;
}
close IN;

my $hutongbiao="/PUBLIC/pipline/database/sheet_feature/";
my ($file,$file2) = judge_lastest_sitedb($hutongbiao);
`perl /PUBLIC/pipline/script/Report/XKFW/V2_171130/script/spectrum_XLSX.feature.pl $file` unless(-e $file2);

my $n=0;
my %index;
my ($line1,$line2,$line2_drug,$line3,$line4,$line5,$line_all)=('','','','','','');
my $qita_GXY='';
my $qita_CS ='';
my $qita_DM ='';
my $qita_OGTT ='';
my $qita_bingfa='';

open IN, $file2 or die $!;
while(<IN>)
{
	chomp;
	@cut=split /\t/, $_;
	($line1,$line2,$line2_drug,$line3,$line4,$line5,$line_all)=('','','','','','','');
	($qita_GXY,$qita_CS,$qita_DM,$qita_OGTT,$qita_bingfa)=('','','','','');
		
	if($n==0)
	{
		for(my $i=0;$i<@cut;$i++)
		{
			$index{$cut[$i]}=$i;
		}
	}
	else{
		next if (!exists $id_hash{$cut[ $index{"样本编号"} ]});
		
		if($version eq "v1" && exists $cut[ $index{"全外怀疑/关心疾病+临床诊断及相关症状"} ] && $cut[ $index{"全外怀疑/关心疾病+临床诊断及相关症状"} ] ne "-" )
		{
			$line1.="疾病及临床诊断: ".$cut[ $index{"全外怀疑/关心疾病+临床诊断及相关症状"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"临床诊断及备注"} ] && $cut[ $index{"临床诊断及备注"} ] ne "-")
		{
			$line1.="临床诊断及备注: ".$cut[ $index{"临床诊断及备注"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"全外怀疑/关心疾病（1-2种）"} ] &&  $cut[ $index{"全外怀疑/关心疾病（1-2种）"} ] ne "-" ) 
		{
			$line1.="全外怀疑/关心疾病: ".$cut[ $index{"全外怀疑/关心疾病（1-2种）"} ].", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"临床诊断疾病"} ] && $cut[ $index{"临床诊断疾病"} ] ne "-" ) 
		{
			$line1.="临床诊断疾病: ".$cut[ $index{"临床诊断疾病"} ].", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"临床症状及发病年龄"} ] && $cut[ $index{"临床症状及发病年龄"} ] ne "-" )
		{
			$line1.="临床症状及发病年龄: ".$cut[ $index{"临床症状及发病年龄"} ]."; \\ \\ ";
		}
		
		##### 高血压 #####
		if(exists $cut[ $index{"高血压性质"} ] && $cut[ $index{"高血压性质"} ] ne "-")
		{
			$line2.="高血压性质: ".$cut[ $index{"高血压性质"} ].", \\ \\ ";
		}
		if(exists $cut[ $index{"发病时间-年份"} ] && $cut[ $index{"发病时间-年份"} ] ne "-")
		{
			$line2.="发病时间: ".$cut[ $index{"发病时间-年份"} ].", \\ \\ ";
		}
		if(exists $cut[ $index{"首发年龄"} ] && $cut[ $index{"首发年龄"} ] ne "-")
		{
			$line2.="首发年龄: ".$cut[ $index{"首发年龄"} ].", \\ \\ ";
		}
		if(exists $cut[ $index{"病程中血压最高达 （mmHg）"} ] && $cut[ $index{"病程中血压最高达 （mmHg）"} ] ne "-")
		{
			$line2.="病程中血压最高达 （mmHg）: ".$cut[ $index{"病程中血压最高达 （mmHg）"} ].", \\ \\ ";
		}

		if(exists $cut[ $index{"目前服用几种药物"} ] && $cut[ $index{"目前服用几种药物"} ] ne "-")
		{
			$line2_drug.="目前服用药物: ".$cut[ $index{"目前服用几种药物"} ]."种, \\ \\ ";
		}
		if($cut[ $index{"CCB"} ] ne "-" )
		{
			$line2_drug.="CCB,";
		}
		if($cut[ $index{"ACEI"} ] ne "-" )
		{
			$line2_drug.="ACEI,";
		}
		if($cut[ $index{"ARB"} ] ne "-" )
		{
			$line2_drug.="ARB,";
		}
		if($cut[ $index{"β-blocker"} ] ne "-" )
		{
			$line2_drug.="β-blocker,";
		}
		if($cut[ $index{"α-blocker"} ] ne "-" )
		{
			$line2_drug.="α-blocker,";
		}
		if($cut[ $index{"氢氯噻嗪"} ] ne "-" )
		{
			$line2_drug.="氢氯噻嗪,";
		}
		if($cut[ $index{"螺内酯"} ] ne "-" )
		{
			$line2_drug.="螺内酯,";
		}
		if($cut[ $index{"阿米洛利"} ] ne "-" )
		{
			$line2_drug.="阿米洛利,";
		}
		if($cut[ $index{"其他"} ] ne "-" )
		{
			$line2_drug.=$cut[ $index{"其他"} ]."; ";
		}
		
		if($line2_drug ne '' && $line2_drug=~/目前服用药物/)
		{
			$line2.=$line2_drug." \\ \\";
		}
		elsif($line2_drug ne '' && $line2_drug!~/目前服用药物/)
		{
			$line2.="目前服用药物：".$line2_drug." \\ \\";
		}

		if($cut[ $index{"心悸、头痛、多汗"} ] ne "-" )
		{
			$qita_GXY.="心悸、头痛、多汗," if($cut[ $index{"心悸、头痛、多汗"} ] eq "√");
			$qita_GXY.=$cut[ $index{"心悸、头痛、多汗"} ].", " if($cut[ $index{"心悸、头痛、多汗"} ] ne "√");
		}
		if($cut[ $index{"多发性腹痛"} ] ne "-" )
		{
			$qita_GXY.="多发性腹痛,";
		}
		if($cut[ $index{"代谢性碱中毒"} ] ne "-" )
		{
			$qita_GXY.="代谢性碱中毒,";
		}
		if($cut[ $index{"癫痫"} ] ne "-" )
		{
			$qita_GXY.="癫痫,";
		}
		if($cut[ $index{"牛奶咖啡斑"} ] ne "-" )
		{
			$qita_GXY.="牛奶咖啡斑,";
		}
		if($cut[ $index{"假两性畸形"} ] ne "-" )
		{
			$qita_GXY.="假两性畸形,";
		}
		if($version eq "v2" && exists $cut[ $index{"其他症状"} ] && $cut[ $index{"其他症状"} ] ne "-" )
		{
			$qita_GXY.="其他:".$cut[ $index{"其他症状"} ]."; \\ \\ ";
		}
		$line2.="其他症状: ".$qita_GXY." \\ \\ "  if($qita_GXY ne '' );
		
		if($cut[ $index{"是否合并低血钾"} ] ne "-" )
		{
			$line2.="是否合并低血钾: ".$cut[ $index{"是否合并低血钾"} ].", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"有无补钾"} ] && $cut[ $index{"有无补钾"} ] ne "-" )
		{
			$line2.="有无补钾: ".$cut[ $index{"有无补钾"} ].", \\ \\ ";
		}
		if($cut[ $index{"平时血钾（mmol/L）"} ] ne "-" )
		{
			$line2.="平时血钾（mmol/L）: ".$cut[ $index{"平时血钾（mmol/L）"} ].", \\ \\ ";
		}
		if($cut[ $index{"病程中最低血钾"} ] ne "-" )
		{
			$line2.="病程中最低血钾: ".$cut[ $index{"病程中最低血钾"} ]."; \\ \\ ";
		}
		if($cut[ $index{"肾上腺CT"} ] ne "-" )
		{
			$line2.="肾上腺CT: ".$cut[ $index{"肾上腺CT"} ].", \\ \\ ";
		}
		if($cut[ $index{"具体描述"} ] ne "-" )
		{
			$line2.="肾上腺CT-具体描述: ".$cut[ $index{"具体描述"} ]."; \\ \\ ";
		}
		
		if($version eq "v1" && exists $cut[ $index{"体位（卧位/立位）"} ] && $cut[ $index{"体位（卧位/立位）"} ] ne "-" )
		{
			$line2.=$cut[ $index{"体位（卧位/立位）"} ].", \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"肾素活性 ng/ml/h  高/低？"} ] && $cut[ $index{"肾素活性 ng/ml/h  高/低？"} ] ne "-" )
		{
			$line2.="肾素活性 ng/ml/h: ".$cut[ $index{"肾素活性 ng/ml/h  高/低？"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"直接肾素 ng/ml  高/低？"} ] && $cut[ $index{"直接肾素 ng/ml  高/低？"} ] ne "-" )
		{
			$line2.="直接肾素 ng/ml: ".$cut[ $index{"直接肾素 ng/ml  高/低？"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"醛固酮 ng/dl  高/低？"} ] && $cut[ $index{"醛固酮 ng/dl  高/低？"} ] ne "-" )
		{
			$line2.="醛固酮 ng/dl: ".$cut[ $index{"醛固酮 ng/dl  高/低？"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"血皮质醇（ug/dl）"} ] && $cut[ $index{"血皮质醇（ug/dl）"} ] ne "-" )
		{
			$line2.="血皮质醇（ug/dl）: ".$cut[ $index{"血皮质醇（ug/dl）"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"ACTH（pg/ml）"} ] && $cut[ $index{"ACTH（pg/ml）"} ] ne "-" )
		{
			$line2.="ACTH（pg/ml）: ".$cut[ $index{"ACTH（pg/ml）"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"血脂水平（正常/异常？）"} ] && $cut[ $index{"血脂水平（正常/异常？）"} ] ne "-" )
		{
			$line2.="血脂水平: ".$cut[ $index{"血脂水平（正常/异常？）"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"血/尿儿茶酚胺（正常/异常？）"} ] && $cut[ $index{"血/尿儿茶酚胺（正常/异常？）"} ] ne "-" )
		{
			$line2.="血/尿儿茶酚胺: ".$cut[ $index{"血/尿儿茶酚胺（正常/异常？）"} ]."; \\ \\ ";
		}
		
		if($version eq "v2" && exists $cut[ $index{"卧位-肾素活性-检测结果 ng/ml/h？  高/低？（参考范围）"} ] && $cut[ $index{"卧位-肾素活性-检测结果 ng/ml/h？  高/低？（参考范围）"} ] ne "-" )
		{
			$line2.="卧位-肾素活性 ng/ml/h: ".&sub_float($cut[ $index{"卧位-肾素活性-检测结果 ng/ml/h？  高/低？（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"立位-肾素活性-检测结果 ng/ml/h？  高/低？（参考范围）"} ] && $cut[ $index{"立位-肾素活性-检测结果 ng/ml/h？  高/低？（参考范围）"} ] ne "-" )
		{
			$line2.="立位-肾素活性 ng/ml/h: ".&sub_float($cut[ $index{"立位-肾素活性-检测结果 ng/ml/h？  高/低？（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"卧位-醛固酮pg/ml  高/低？（参考范围）"} ] && $cut[ $index{"卧位-醛固酮pg/ml  高/低？（参考范围）"} ] ne "-" )
		{
			$line2.="卧位-醛固酮pg/ml: ".&sub_float($cut[ $index{"卧位-醛固酮pg/ml  高/低？（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"立位-醛固酮pg/ml  高/低？"} ] && $cut[ $index{"立位-醛固酮pg/ml  高/低？"} ] ne "-" )
		{
			$line2.="立位-醛固酮pg/ml: ".&sub_float($cut[ $index{"立位-醛固酮pg/ml  高/低？"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"卧位-直接肾素 ng/ml  高/低？（参考范围）"} ] && $cut[ $index{"卧位-直接肾素 ng/ml  高/低？（参考范围）"} ] ne "-" )
		{
			$line2.="卧位-直接肾素 ng/ml: ".&sub_float($cut[ $index{"卧位-直接肾素 ng/ml  高/低？（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"立位-直接肾素 ng/ml  高/低？（参考范围）"} ] && $cut[ $index{"立位-直接肾素 ng/ml  高/低？（参考范围）"} ] ne "-" )
		{
			$line2.="立位-直接肾素 ng/ml: ".&sub_float($cut[ $index{"立位-直接肾素 ng/ml  高/低？（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血皮质醇（nmol/L）-0:00 （参考范围）"} ] && $cut[ $index{"血皮质醇（nmol/L）-0:00 （参考范围）"} ] ne "-" )
		{
			$line2.="血皮质醇（nmol/L）-0:00: ".&sub_float($cut[ $index{"血皮质醇（nmol/L）-0:00 （参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血皮质醇（nmol/L）-8:00（参考范围）"} ] && $cut[ $index{"血皮质醇（nmol/L）-8:00（参考范围）"} ] ne "-" )
		{
			$line2.="血皮质醇（nmol/L）-8:00: ".&sub_float($cut[ $index{"血皮质醇（nmol/L）-8:00（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血皮质醇（nmol/L）-16:00（参考范围）"} ] && $cut[ $index{"血皮质醇（nmol/L）-16:00（参考范围）"} ] ne "-" )
		{
			$line2.="血皮质醇（nmol/L）-16:00: ".&sub_float($cut[ $index{"血皮质醇（nmol/L）-16:00（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"ACTH（pg/ml）-0:00（参考范围）"} ] && $cut[ $index{"ACTH（pg/ml）-0:00（参考范围）"} ] ne "-" )
		{
			$line2.="ACTH（pg/ml）-0:00: ".&sub_float($cut[ $index{"ACTH（pg/ml）-0:00（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"ACTH（pg/ml）-8:00（参考范围）"} ] && $cut[ $index{"ACTH（pg/ml）-8:00（参考范围）"} ] ne "-" )
		{
			$line2.="ACTH（pg/ml）-8:00: ".&sub_float($cut[ $index{"ACTH（pg/ml）-8:00（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"ACTH（pg/ml）-16:00-参考范围"} ] && $cut[ $index{"ACTH（pg/ml）-16:00-参考范围"} ] ne "-" )
		{
			$line2.="ACTH（pg/ml）-16:00: ".&sub_float($cut[ $index{"ACTH（pg/ml）-16:00-参考范围"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血儿茶酚胺（正常/异常？）ug/dl （参考范围）"} ] && $cut[ $index{"血儿茶酚胺（正常/异常？）ug/dl （参考范围）"} ] ne "-" )
		{
			$line2.="血儿茶酚胺 ug/dl: ".&sub_float($cut[ $index{"血儿茶酚胺（正常/异常？）ug/dl （参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"尿儿茶酚胺（正常/异常？）ug/dl（参考范围）"} ] && $cut[ $index{"尿儿茶酚胺（正常/异常？）ug/dl（参考范围）"} ] ne "-" )
		{
			$line2.="尿儿茶酚胺 ug/dl: ".&sub_float($cut[ $index{"尿儿茶酚胺（正常/异常？）ug/dl（参考范围）"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血钠mmol/L（参考范围）"} ] && $cut[ $index{"血钠mmol/L（参考范围）"} ] ne "-" )
		{
			$line2.="血钠mmol/L: ".$cut[ $index{"血钠mmol/L（参考范围）"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血镁mmol/l（参考范围）"} ] && $cut[ $index{"血镁mmol/l（参考范围）"} ] ne "-" )
		{
			$line2.="血镁mmol/l: ".$cut[ $index{"血镁mmol/l（参考范围）"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血氯mmol/l（参考范围）"} ] && $cut[ $index{"血氯mmol/l（参考范围）"} ] ne "-" )
		{
			$line2.="血氯mmol/l: ".$cut[ $index{"血氯mmol/l（参考范围）"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血钙mmol/l（参考范围）"} ] && $cut[ $index{"血钙mmol/l（参考范围）"} ] ne "-" )
		{
			$line2.="血钙mmol/l: ".$cut[ $index{"血钙mmol/l（参考范围）"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血脂水平（正常/异常？）"} ] && $cut[ $index{"血脂水平（正常/异常？）"} ] ne "-" )
		{
			$line2.="血脂水平: ".$cut[ $index{"血脂水平（正常/异常？）"} ]."; \\ \\ ";
		}	
		
		if($cut[ $index{"性激素（正常/异常？）"} ] ne "-" )
		{
			$line2.="性激素: ".$cut[ $index{"性激素（正常/异常？）"} ]."; \\ \\ ";
		}
		if($cut[ $index{"其他激素或生化指标：指标名称及水平"} ] ne "-" )
		{
			$line2.="其他激素或生化指标: ".$cut[ $index{"其他激素或生化指标：指标名称及水平"} ]."; \\ \\ ";
		}
		
		##### 猝死 #####
		if($cut[ $index{"心悸"} ] ne "-" )
		{
			$qita_CS.="心悸,";
		}
		if($cut[ $index{"气短"} ] ne "-" )
		{
			$qita_CS.="气短,";
		}
		if($cut[ $index{"运动后气短"} ] ne "-" )
		{
			$qita_CS.="运动后气短,";
		}
		if($cut[ $index{"胸闷"} ] ne "-" )
		{
			$qita_CS.="胸闷,";
		}
		if($cut[ $index{"晕厥"} ] ne "-" )
		{
			$qita_CS.="晕厥,";
		}
		if($cut[ $index{"心衰"} ] ne "-" )
		{
			$qita_CS.="心衰,";
		}
		if($cut[ $index{"冠心病"} ] ne "-" )
		{
			$qita_CS.="冠心病,";
		}
		if($cut[ $index{"黄色瘤"} ] ne "-" )
		{
			$qita_CS.="黄色瘤,";
		}
		if($cut[ $index{"黄踺瘤"} ] ne "-" )
		{
			$qita_CS.="黄踺瘤,";
		}
		if($cut[ $index{"心梗"} ] ne "-" )
		{
			$qita_CS.="心梗,";
		}
		if($cut[ $index{"脑卒中"} ] ne "-" )
		{
			$qita_CS.="脑卒中,";
		}
		if($cut[ $index{"房颤"} ] ne "-" )
		{
			$qita_CS.="房颤,";
		}
		$line3.="临床症状: ".$qita_CS." \\ \\ "  if($qita_CS ne '' );
		
		
		if($cut[ $index{"心脏B超"} ] ne "-" )
		{
			$line3.="心脏B超: ".$cut[ $index{"心脏B超"} ]."; \\ \\ ";
		}
		if($cut[ $index{"彩超"} ] ne "-" )
		{
			$line3.="彩超: ".$cut[ $index{"彩超"} ]."; \\ \\ ";
		}
		if($cut[ $index{"心动图"} ] ne "-" )
		{
			$line3.="心动图: ".$cut[ $index{"心动图"} ]."; \\ \\ ";
		}
		if($cut[ $index{"心电图"} ] ne "-" )
		{
			$line3.="心电图: ".$cut[ $index{"心电图"} ]."; \\ \\ ";
		}
		if($cut[ $index{"冠脉造影"} ] ne "-" )
		{
			$line3.="冠脉造影: ".$cut[ $index{"冠脉造影"} ]."; \\ \\ ";
		}
		if($cut[ $index{"心脏MRI"} ] ne "-" )
		{
			$line3.="心脏MRI: ".$cut[ $index{"心脏MRI"} ]."; \\ \\ ";
		}
		if($cut[ $index{"血脂水平"} ] ne "-" )
		{
			$line3.="血脂水平: ".$cut[ $index{"血脂水平"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"LDL-C水平"} ] && $cut[ $index{"LDL-C水平"} ] ne "-" )
		{
			$line3.="LDL-C水平: ".&sub_float($cut[ $index{"LDL-C水平"} ])."; \\ \\ ";
		}
		
		##### 糖尿病 #####
		if($version eq "v2" && exists $cut[ $index{"发病年龄"} ] && $cut[ $index{"发病年龄"} ] ne "-" )
		{
			$line4.="发病年龄: ".$cut[ $index{"发病年龄"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"BMI指数"} ] && $cut[ $index{"BMI指数"} ] ne "-" )
		{
			$line4.="BMI指数: ".$cut[ $index{"BMI指数"} ]."; \\ \\ ";
		}
		if($cut[ $index{"多饮"} ] ne "-" )
		{
			$qita_DM.="多饮,";
		}
		if($cut[ $index{"多尿"} ] ne "-" )
		{
			$qita_DM.="多尿,";
		}
		if($cut[ $index{"多食"} ] ne "-" )
		{
			$qita_DM.="多食,";
		}
		if($cut[ $index{"体重减轻"} ] ne "-" )
		{
			$qita_DM.="体重减轻,";
		}
		if($cut[ $index{"酮症酸中毒"} ] ne "-" )
		{
			$qita_DM.="酮症酸中毒,";
		}
		if($cut[ $index{"乳酸酸中毒"} ] ne "-" )
		{
			$qita_DM.="乳酸酸中毒,";
		}
		if($cut[ $index{"低血糖"} ] ne "-" )
		{
			$qita_DM.="低血糖,";
		}
		if($cut[ $index{"肥胖"} ] ne "-" )
		{
			$qita_DM.="肥胖,";
		}
		if($cut[ $index{"肾脏异常表现"} ] ne "-" )
		{
			$qita_DM.="肾脏异常表现,";
		}
		if($cut[ $index{"生殖器异常"} ] ne "-" )
		{
			$qita_DM.="生殖器异常,";
		}
		if($version eq "v2" && exists $cut[ $index{"妊娠糖尿病史"} ] && $cut[ $index{"妊娠糖尿病史"} ] ne "-" )
		{
			$qita_DM.="妊娠糖尿病史,";
		}
		if($version eq "v2" && exists $cut[ $index{"脑病"} ] && $cut[ $index{"脑病"} ] ne "-" )
		{
			$qita_DM.="脑病,";
		}
		if($version eq "v2" && exists $cut[ $index{"耳聋"} ] && $cut[ $index{"耳聋"} ] ne "-" )
		{
			$qita_DM.="耳聋,";
		}
		if($version eq "v2" && exists $cut[ $index{"心肌病"} ] && $cut[ $index{"心肌病"} ] ne "-" )
		{
			$qita_DM.="心肌病,";
		}
		if($version eq "v2" && exists $cut[ $index{"黑棘皮病"} ] && $cut[ $index{"黑棘皮病"} ] ne "-" )
		{
			$qita_DM.="黑棘皮病,";
		}
		$line4.="临床症状: ".$qita_DM." \\ \\ "  if($qita_DM ne '' );
		
		if($cut[ $index{"胰岛免疫抗体（IAA、IA2、GAD）"} ] ne "-" )
		{
			$line4.="胰岛免疫抗体（IAA、IA2、GAD）: ".$cut[ $index{"胰岛免疫抗体（IAA、IA2、GAD）"} ]."; \\ \\ ";
		}
		if($cut[ $index{"DM血脂水平"} ] ne "-" )
		{
			$line4.="血脂水平: ".$cut[ $index{"DM血脂水平"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"C-肽"} ] && $cut[ $index{"C-肽"} ] ne "-" )
		{
			$line4.="C-肽: ".$cut[ $index{"C-肽"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"糖化血红蛋白（HbAIc）"} ] && $cut[ $index{"糖化血红蛋白（HbAIc）"} ] ne "-" )
		{
			$line4.="糖化血红蛋白（HbAIc）: ".$cut[ $index{"糖化血红蛋白（HbAIc）"} ]."; \\ \\ ";
		}
		if($version eq "v1" && exists $cut[ $index{"空腹血糖"} ] && $cut[ $index{"空腹血糖"} ] ne "-" )
		{
			$line4.="空腹血糖: ".$cut[ $index{"空腹血糖"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C-肽-是否服药"} ] && $cut[ $index{"C-肽-是否服药"} ] ne "-" )
		{
			$line4.="C-肽: ".$cut[ $index{"C-肽-是否服药"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C肽-空腹"} ] && $cut[ $index{"C肽-空腹"} ] ne "-" )
		{
			$line4.="C肽-空腹: ".&sub_float($cut[ $index{"C肽-空腹"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C肽-餐后2h"} ] && $cut[ $index{"C肽-餐后2h"} ] ne "-" )
		{
			$line4.="C肽-餐后2h: ".&sub_float($cut[ $index{"C肽-餐后2h"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"糖化血红蛋白（HbA1c）-是否服药"} ] && $cut[ $index{"糖化血红蛋白（HbA1c）-是否服药"} ] ne "-" )
		{
			$line4.="糖化血红蛋白（HbA1c）: ".$cut[ $index{"糖化血红蛋白（HbA1c）-是否服药"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"糖化血红蛋白（HbA1c）-空腹"} ] && $cut[ $index{"糖化血红蛋白（HbA1c）-空腹"} ] ne "-" )
		{
			$line4.="糖化血红蛋白（HbA1c）-空腹: ".$cut[ $index{"糖化血红蛋白（HbA1c）-空腹"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"糖化血红蛋白（HbA1c）-餐后2h"} ] && $cut[ $index{"糖化血红蛋白（HbA1c）-餐后2h"} ] ne "-" )
		{
			$line4.="糖化血红蛋白（HbA1c）-餐后2h: ".$cut[ $index{"糖化血红蛋白（HbA1c）-餐后2h"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-是否服药"} ] && $cut[ $index{"血糖-是否服药"} ] ne "-" )
		{
			$line4.="血糖: ".$cut[ $index{"血糖-是否服药"} ]."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-空腹"} ] && $cut[ $index{"血糖-空腹"} ] ne "-" )
		{
			$line4.="血糖-空腹: ".&sub_float($cut[ $index{"血糖-空腹"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-餐后2h"} ] && $cut[ $index{"血糖-餐后2h"} ] ne "-" )
		{
			$line4.="血糖-餐后2h: ".&sub_float($cut[ $index{"血糖-餐后2h"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C肽-0h"} ] && $cut[ $index{"C肽-0h"} ] ne "-" )
		{
			$qita_OGTT.="C肽-0h: ".&sub_float($cut[ $index{"C肽-0h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C肽-0.5h"} ] && $cut[ $index{"C肽-0.5h"} ] ne "-" )
		{
			$qita_OGTT.="C肽-0.5h: ".&sub_float($cut[ $index{"C肽-0.5h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C肽-1h"} ] && $cut[ $index{"C肽-1h"} ] ne "-" )
		{
			$qita_OGTT.="C肽-1h: ".&sub_float($cut[ $index{"C肽-1h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C肽-2h"} ] && $cut[ $index{"C肽-2h"} ] ne "-" )
		{
			$qita_OGTT.="C肽-2h: ".&sub_float($cut[ $index{"C肽-2h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"C肽-3h"} ] && $cut[ $index{"C肽-3h"} ] ne "-" )
		{
			$qita_OGTT.="C肽-3h: ".&sub_float($cut[ $index{"C肽-3h"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"胰岛素-0h"} ] && $cut[ $index{"胰岛素-0h"} ] ne "-" )
		{
			$qita_OGTT.="胰岛素-0h: ".&sub_float($cut[ $index{"胰岛素-0h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"胰岛素-0.5h"} ] && $cut[ $index{"胰岛素-0.5h"} ] ne "-" )
		{
			$qita_OGTT.="胰岛素-0.5h: ".&sub_float($cut[ $index{"胰岛素-0.5h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"胰岛素-1h"} ] && $cut[ $index{"胰岛素-1h"} ] ne "-" )
		{
			$qita_OGTT.="胰岛素-1h: ".&sub_float($cut[ $index{"胰岛素-1h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"胰岛素-2h"} ] && $cut[ $index{"胰岛素-2h"} ] ne "-" )
		{
			$qita_OGTT.="胰岛素-2h: ".&sub_float($cut[ $index{"胰岛素-2h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"胰岛素-3h"} ] && $cut[ $index{"胰岛素-3h"} ] ne "-" )
		{
			$qita_OGTT.="胰岛素-3h: ".&sub_float($cut[ $index{"胰岛素-3h"} ])."; \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-0h"} ] && $cut[ $index{"血糖-0h"} ] ne "-" )
		{
			$qita_OGTT.="血糖-0h: ".&sub_float($cut[ $index{"血糖-0h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-0.5h"} ] && $cut[ $index{"血糖-0.5h"} ] ne "-" )
		{
			$qita_OGTT.="血糖-0.5h: ".&sub_float($cut[ $index{"血糖-0.5h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-1h"} ] && $cut[ $index{"血糖-1h"} ] ne "-" )
		{
			$qita_OGTT.="血糖-1h: ".&sub_float($cut[ $index{"血糖-1h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-2h"} ] && $cut[ $index{"血糖-2h"} ] ne "-" )
		{
			$qita_OGTT.="血糖-2h: ".&sub_float($cut[ $index{"血糖-2h"} ]).", \\ \\ ";
		}
		if($version eq "v2" && exists $cut[ $index{"血糖-3h"} ] && $cut[ $index{"血糖-3h"} ] ne "-" )
		{
			$qita_OGTT.="血糖-3h: ".&sub_float($cut[ $index{"血糖-3h"} ])."; \\ \\ ";
		}
		$line4.="OGTT实验: ".$qita_OGTT  if($qita_OGTT ne '' );
		
		if($version eq "v2" && exists $cut[ $index{"使用药物"} ] && $cut[ $index{"使用药物"} ] ne "-" )
		{
			$line4.="使用药物: ".$cut[ $index{"使用药物"} ]."; \\ \\ ";
		}
		if($cut[ $index{"是否出现并发症"} ] ne "-" )
		{
			$line4.="是否出现并发症: ".$cut[ $index{"是否出现并发症"} ].", \\ \\ ";
		}
		if($cut[ $index{"糖尿病肾病"} ] ne "-" )
		{
			$qita_bingfa.="糖尿病肾病,";
		}
		if($cut[ $index{"视网膜病变"} ] ne "-" )
		{
			$qita_bingfa.="视网膜病变,";
		}
		if($cut[ $index{"糖尿病病足"} ] ne "-" )
		{
			$qita_bingfa.="糖尿病病足,";
		}
		if($cut[ $index{"其他并发症"} ] ne "-" )
		{
			$qita_bingfa.="其他: ".$cut[ $index{"其他并发症"} ]."; \\ \\ ";
		}
		$line4.=$qita_bingfa." \\ \\ "  if($qita_bingfa ne '' );
		
		##### 其他病史 #####
		if($cut[ $index{"临床诊断"} ] ne "-" )
		{
			$line5.="临床诊断: ".$cut[ $index{"临床诊断"} ].", \\ \\ ";
		}
		if($cut[ $index{"家族史-（高血压/心脏病/糖尿病-----有/无）"} ] ne "-" )
		{
			$line5.="家族史（高血压/心脏病/糖尿病）: ".$cut[ $index{"家族史-（高血压/心脏病/糖尿病-----有/无）"} ].", \\ \\ ";
		}
		if($cut[ $index{"何人（与患者关系）"} ] ne "-" )
		{
			$line5.="何人（与患者关系）: ".$cut[ $index{"何人（与患者关系）"} ].", \\ \\ ";
		}
		if($cut[ $index{"起病年龄"} ] ne "-" )
		{
			$line5.="起病年龄: ".$cut[ $index{"起病年龄"} ].", \\ \\ ";
		}
		if($cut[ $index{"病程"} ] ne "-" )
		{
			$line5.="病程: ".$cut[ $index{"病程"} ].", \\ \\ ";
		}
		if($cut[ $index{"疾病名称"} ] ne "-" )
		{
			$line5.="疾病名称: ".$cut[ $index{"疾病名称"} ].", \\ \\ ";
		}
		if($cut[ $index{"血压"} ] ne "-" )
		{
			$line5.="血压: ".$cut[ $index{"血压"} ].", \\ \\ ";
		}
		if($cut[ $index{"血钾"} ] ne "-" )
		{
			$line5.="血钾: ".$cut[ $index{"血钾"} ].", \\ \\ ";
		}
		if($cut[ $index{"血糖"} ] ne "-" )
		{
			$line5.="血糖: ".$cut[ $index{"血糖"} ].", \\ \\ ";
		}
		if($cut[ $index{"是否有其他家族遗传病史"} ] ne "-" )
		{
			$line5.="其他家族遗传病史: ".$cut[ $index{"是否有其他家族遗传病史"} ]."; \\ \\ ";
		}
		if($cut[ $index{"家族中何人患有何病"} ] ne "-" )
		{
			$line5.="家族中何人患有何病: ".$cut[ $index{"家族中何人患有何病"} ]."; \\ \\ ";
		}

		$line_all=$cut[ $index{"样本编号"} ]."\t";		
		if($line1 ne ''){$line1=~s/,$//;$line1=~s/, \\ \\ $//;$line1=~s/; \\ \\ $//; $line_all.= $line1." \\\\ ";}
		if($line2 ne ''){$line2=~s/,$//;$line2=~s/, \\ \\ $//;$line2=~s/; \\ \\ $//; $line_all.=  "{\\sym{高血压:}} ".$line2." \\\\ ";}
		if($line3 ne ''){$line3=~s/,$//;$line3=~s/, \\ \\ $//;$line3=~s/; \\ \\ $//; $line_all.=  "{\\sym{遗传性心血管病:}} ".$line3." \\\\ ";}
		if($line4 ne ''){$line4=~s/,$//;$line4=~s/, \\ \\ $//;$line4=~s/; \\ \\ $//; $line_all.=  "{\\sym{糖尿病:}} ".$line4." \\\\ ";}
		if($line5 ne ''){$line5=~s/,$//;$line5=~s/, \\ \\ $//;$line5=~s/; \\ \\ $//; $line_all.=  "{\\sym{其他病史:}} ".$line5." ";}
		$line_all=~s/\\\\ $//;
		$line_all=~s/\%/\\\%/g;
		$line_all=~s/\_/\\\_/g;
		print $line_all."\n";
	}
	$n++;

	
}
close IN;

sub sub_float{
	my $input=shift;
	my $out=$input;

	if($input=~/^[-+]?[0-9]*\.?([0-9]+)$/)
	{
		if( (length $1)>3 )
		{
			$out = sprintf "%.3f" , $input; 	
		}
	}

	return $out;
}

sub judge_lastest_sitedb{
        my $product_site=shift;
        my @all_date;
        my $lastest_date;
        my $date;

        my @all_path=glob "$product_site/*";
        foreach(@all_path)
        {
                $date=(split /\//,$_)[-1];
                if($date=~/心康-临床样本信息单-(\d+)\.xlsx/)
                {
                        push @all_date,$1;
                }
        }

        $lastest_date= max(@all_date);
        my $lastest_sitedb_path ="$product_site/心康-临床样本信息单-$lastest_date.xlsx";
        my $lastest_sitedb_path2="$product_site/心康-临床样本信息单-$lastest_date.txt";
        print STDERR "$lastest_sitedb_path\n";
        return $lastest_sitedb_path,$lastest_sitedb_path2;
}
