indir=$1
proj=${indir:5:6}
year="20"${indir:12:2}
month=${indir:14:2}
if [ -z $indir ];then
echo "sh work.sh indir"
exit
fi
if [ "$proj" == "proj06" ];then
	year="01"
	month="全外基础版位点数据"
        ln -s $indir/CNV/result/*/*.merge.cnv.xlsx $indir/report/check/
fi
if [ `cat /XK1/fast_start/Monitor_xinkang/sample_path.file |grep -c $indir` -lt 1 ];then
	perl /XK1/$proj/sample.list/step1.sample.merge.pl /XK1/$proj/
	python /XK1/fast_start/Report_xinkang/merge_sample_path.py
fi
outdir=`cat /XK1/fast_start/Monitor_xinkang/sample_path.file |grep  $indir|cut -f1`
rm $indir/report/check/samples_info
for sample in `cut -f1 $indir/report/01products/cfg.txt |grep -v 样本编号|sort`
do
	for sam in `awk '$2=="'$sample'" {print $3}' /PUBLIC/pipline/database/sheet_hutong/2019心康送样信息单-1.family.txt`
	do 
		awk '$3=="'$sam'" {print $0}' /PUBLIC/pipline/database/sheet_hutong/2019心康送样信息单-1.txt >>$indir/report/check/samples_info
	done		
	#echo $sample
	#awk '$2=="'$sample'" {print $3"\t"$4"\t"$1}' /PUBLIC/pipline/database/sheet_hutong/2019心康送样信息单-1.family.txt >>$indir/report/check/samples
done
echo """source /volume1/homes/wangpenghui/.bashrc
if [ ! -d /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month ];then
mkdir /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month
fi
mkdir /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month/$outdir
mkdir /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month/$outdir/check_newpip
mkdir /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month/$outdir/HgmdClivar_site
cpr $indir/Stat_QC.xls /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month/$outdir
cpr $indir/report/check/* /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month/$outdir/check_newpip
cpr $indir/report/check/HgmdClivar_site/* /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month/$outdir/HgmdClivar_site""" >$indir/report/cpr.sh
echo "python /volume1/production-operations/生产-运营/.log/manageJieduTime.py update -p /volume1/production/03.产品及报告/01.测序产品/02.下机数据/$year.$month/$outdir/check_newpip/samples_info" >>$indir/report/cpr.sh
echo "$indir/report/cpr.sh"
