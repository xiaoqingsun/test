path=$1
cfg=$path/report/01products/cfg.txt
#echo $path,$cfg
for sample in `cut -f1 $cfg|grep -v "样本编号"`
do 
    head -1 $path/cr/$sample/$sample/02.Aln/Stat/$sample\.target.low.depth.gene.bed.info >$path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info
    for gene in `less $path/report/LOG/$sample\.genelist`
        do 
            awk -F'\t' '$9=="'$gene'" {print $0}' $path/cr/$sample/$sample/02.Aln/Stat/$sample\.target.low.depth.gene.bed.info >>$path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info
        done
    python2 /PUBLIC/home/wangpenghui/Findgene/txt2xls.py $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.info
    ln -s $path/cr/$sample/$sample/02.Aln/Stat/$sample\.filter.target.low.depth.gene.bed.xls $path/report/check/
    echo $sample
done
