path=$1
cfg=$path/report/01products/cfg.txt
productsam=`cut -f12 $cfg|grep -v "样本类型"|head -1`

#echo $path,$cfg
if [ `ls $path/report/LOG/ |grep -c genelist` -eq 0 ];then
	for sam in `cat $cfg |grep 单基因高血压|cut -f1`
		do 
			cat /PUBLIC/pipline/database/siteFilterDB/product/GXY/GXY_shaicha_gene.list.update >> $path/report/LOG/$sam.genelist
		done
	for sam in `cat $cfg |grep 心源性猝死|cut -f1`
		do 
			cat /PUBLIC/pipline/database/siteFilterDB/product/CuS/cusi_xinkang.geneV4 >> $path/report/LOG/$sam.genelist
		done
	for sam in `cat $cfg |grep 单基因糖尿病|cut -f1`
		do
			cat /PUBLIC/pipline/database/siteFilterDB/product/DM/DM_xinkang.gene >> $path/report/LOG/$sam.genelist
		done
fi
for sample in `cut -f1 $cfg|grep -v "样本编号"`
do 
    head -1 $path/cr/$sample/$sample/02.Aln/Stat/$sample\.target.low.depth.gene.bed.info >$path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info
    for gene in `cat $path/report/LOG/$sample\.genelist|sort|uniq`
        do
            awk -F'\t' '$9=="'$gene'" {print $0}' $path/cr/$sample/$sample/02.Aln/Stat/$sample\.target.low.depth.gene.bed.info >> $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info
        done
    less $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info | grep -v CHR | /PUBLIC/pipline/script/Bin/bedtools/intersectBed -a -  -b /PUBLIC/pipline/database/Clinvar/variant_summary.snp_indel.final.bedtools.txt -loj > $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info.clinvar
    perl /PUBLIC/pipline/script/Report/XKFW/V3_180824/script/filter_low_dp.clinvar_infor.pl $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info.clinvar /PUBLIC/pipline/database/Clinvar/variant_summary.snp.final.txt /PUBLIC/pipline/database/Clinvar/variant_summary.indel.final.txt > $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info.clinvar.filt
    perl /PUBLIC/pipline/script/Report/XKFW/V3_180824/script/supply_info_lowdp.pl  -pos   20,21,22,23  -in $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info.clinvar.filt  -bed $productsam
    python2 /PUBLIC/work/wangpenghui/生产交互表/txt2xlsx.py $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info.clinvar.filt
    ln -s $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info.clinvar.xlsx $path/report/check/
    echo $sample
done
sh /PUBLIC/pipline/script/Report/XKFW/V3_180824/script/release_data.sh ${path:0:23}
