#-*- coding:utf-8 -*-
import os
import sys
reload(sys)
sys.setdefaultencoding("utf-8")
import shutil
import logging
import argparse
import traceback
import time
import decimal
import xlrd
import os
import glob
import sys
import re
from string import Template
from collections import defaultdict


def argInit():
    parser = argparse.ArgumentParser(
        description="report_create program.")
    parser.add_argument(
        '-dir', help='file CS.filter.txt', required=True)
    parser.add_argument(
        '-part', help='2 / 3 / 4', required=True)
    parser.add_argument(
        '-out', help='project name', default=os.getcwd())
    parser.add_argument(
        '-sample', help='sampleid', required=True)
    parser.add_argument(
        '-module', help='module', required=True)
    parser.add_argument(
        '-yinyang', help='阴阳性', required=True)
    argv = vars(parser.parse_args())
    global filterFile, part,outDir,sampleid,module,yinyang_file
    filterFile = argv['dir'].strip()
    part = argv['part'].strip()
    module = argv['module'].strip()
    outDir = os.path.abspath(argv['out'].strip())
    sampleid=argv['sample'].strip()
    yinyang_file=argv['yinyang'].strip()

def get_filter(filterFile,sampleid):
    filter_li=[]
    filter_fu=[]
    filter_bz=[]
    id_li=[]
    if os.path.isfile(filterFile) == False:
        sys.stderr.write('%s is not a file' % filterFile)
        sys.exit(1)
    try:
        with open(filterFile, 'rU') as file_handle:
            for line in file_handle:
                if re.match('^ID',line.strip()):continue
                if re.match('^####',line.strip()):continue
                line_list=[x.strip() for x in line.strip().split('\t')]
	        if re.match('^###',line.strip()):
		    newsam=re.sub('###','',line_list[0])
                    if  newsam==sampleid:
			line_1=re.sub('%','\%',line)
                        line_2=re.sub('_','\_',line_1)
			line_3=re.sub('~','-',line_2) ### panqi 20181213
                        filter_bz.append(line_3.strip())
                elif re.match('^##',line.strip()): 
                    newsam=re.sub('##','',line_list[0])
                    if  newsam==sampleid:
                        line_1=re.sub('%','\%',line)
                        line_2=re.sub('_','\_',line_1)
			line_3=re.sub('~','-',line_2) ### panqi 20181213
                        filter_fu.append(line_3.strip())
                else:
                    id_li.append(line_list[0])
                    if ''.join(line_list[0:5]) == '-----':continue
                    if line_list[0]==sampleid:
                        line_1=re.sub('%','\%',line)
                        line_2=re.sub('_','\_',line_1)
			line_3=re.sub('~','-',line_2) ### panqi 20181213
                        filter_li.append(line_3.strip())
                    else:
                        pass
        if sampleid not in id_li:
            sys.stderr.write('sample %s is not in %s'%(sampleid,filterFile))
            sys.exit(1)
    except Exception as e:
        traceback.print_exc(e)
    return filter_li,filter_fu,filter_bz

def prepare_allDic(files,filesfu,filesbz,fff):
    disease_Dic = defaultdict(list)
    site_Dic = defaultdict(list)
    site_li = []
    disease_li= []
    het_hom ={'het':u'杂合','hom':u'纯合','hem':u'半合子'} #wph add hem 181212
    data_Dic=read_CSdata()
    if fff==1:
        try:
            for line in files:
                line_list=[x.strip() for x in line.strip().split('\t')]
                if len(line_list) < 88:
                    sys.stderr.write('filter文件列数不足，请检查')
                    sys.exit(1)
            ##获取disease_Dic
                gene_dis_key=re.sub(' ','',line_list[85])+":"+ re.sub(' ','',line_list[83]) #gene:disease_subtype作为key(去空格)
                if gene_dis_key not in disease_li:     #排除相同gene（不同位点）有多个相同疾病
                     disease_li.append(gene_dis_key)
                     if data_Dic[gene_dis_key]:
                         detail=data_Dic[gene_dis_key]
                         detail.append(line_list[79])  #加入疾病描述
                         disease_Dic[line_list[85]].append(detail)
                     else:
                         sys.stderr.write(u'请检查filter内%s与CS.data.txt中的gene名、疾病亚型名是否一致'%gene_dis_key)
                         sys.exit(1)
                else:
                    pass
            ##获取site_Dic
                site=line_list[1]+":"+line_list[2]+":"+line_list[3]+":"+line_list[4]+":"+line_list[5] #chr:start:end:ref:alt
                if site not in site_li:    #排除相同基因（不同疾病）有多个相同位点
                    chr_site=get_chr(line_list[1],line_list[2],line_list[3])
                    if line_list[88]=='-':
                        exon='-'
                    else:
	    	        ncbi_tran=line_list[88]   #changed by sunxq
                        exons=ncbi_tran.split(':')
                        exon=exons[0]+'\\\\'+exons[1]
                    mute_type=get_mut_type(line_list[9])  
                    level_value=get_level_value(line_list[78])
                    site_detail=[line_list[86],line_list[87],chr_site,exon,het_hom[line_list[31].lower()],mute_type,line_list[78],line_list[82],line_list[80],line_list[81],level_value]  #nue,acid,chr_site,exon,het/hom/,mute_type,zhibingxing,literature,gene_des,site_des,level_value
                #print line_list[9],line_list[10],chr_site,exon,het_hom[line_list[26].lower()],mute_type,line_list[66],line_list[70],line_list[68],line_list[69]
                    site_Dic[line_list[85]].append(site_detail)
                    site_li.append(site)
                else:
                    pass
        except Exception as e:
                traceback.print_exc(e)
        sort_disease = sort_Dic(disease_Dic)
        sort_site = sort_siteDic(site_Dic)
    else:
        sort_disease = disease_Dic
        sort_site = site_Dic
    disease_Dic_fubiao = defaultdict(list)
    site_Dic_fubiao = defaultdict(list)
    site_li_fubiao = []
    disease_li_fubiao = []
    try:
        for line in filesfu:
            line_list=[x.strip() for x in line.strip().split('\t')]
            genekey=re.sub(' ','',line_list[7])
            disea=re.sub(' ','',line_list[83])
            for gdkey in data_Dic:
                line_list2=[x.strip() for x in gdkey.strip().split(':')]
                if  line_list2[0] == genekey  and disea == line_list2[1]:
                        if gdkey not in disease_li_fubiao:
                            disease_li_fubiao.append(gdkey)
                            detail=data_Dic[gdkey]
                        #detail.append(line_list2[1])
                            disease_Dic_fubiao[genekey].append(detail)
                        #print detail[0]+detail[1]+detail[2]+detail[3]+detail[4]+detail[5]+detail[6]
            site=line_list[1]+":"+line_list[2]+":"+line_list[3]+":"+line_list[4]+":"+line_list[5]
            if site not in site_li_fubiao:
                if line_list[88]=='-':
                    exon='-'
                else:
                    ncbi_tran=line_list[88]   #changed by sunxq
                    exons=ncbi_tran.split(':')
                    exon=exons[0]+'\\\\'+exons[1]
                mute_type=get_mut_type(line_list[9])   
                level_value=get_level_value(line_list[78]) 
                site_detail=[line_list[86],line_list[87],exon,het_hom[line_list[31].lower()],mute_type,line_list[21],line_list[77],line_list[78].replace("临床意义","临床意义\\\\"),level_value]
                site_Dic_fubiao[line_list[7]].append(site_detail)
                site_li_fubiao.append(site)
            else:
                 pass
    except Exception as e:
            traceback.print_exc(e)
    Dic_bz = defaultdict(list)
    try:
	for line in filesbz:
             line_list=[x.strip() for x in line.strip().split('\t')]
             genekey=re.sub(' ','',line_list[7])
             Dic_bz[genekey].append(line_list[83]) 
    except Exception as e:
            traceback.print_exc(e)
    #sort_disease_fubiao = sort_Dic(disease_Dic_fubiao)
    sort_disease_fubiao=disease_Dic_fubiao
    sort_site_fubiao = sort_siteDic_fubiao(site_Dic_fubiao)
    return sort_disease,sort_site,sort_disease_fubiao,sort_site_fubiao,Dic_bz

def sort_siteDic(Dic):
    sorted_Dic={}
    for key in Dic:
        sorted_Dic[key]=sorted(Dic[key],key=lambda x:x[10])
    return sorted_Dic

def sort_siteDic_fubiao(Dic):
    sorted_Dic={}
    for key in Dic:
        sorted_Dic[key]=sorted(Dic[key],key=lambda x:x[8])
    return sorted_Dic

def sort_Dic(Dic):
    # [疾病名，疾病亚型，遗传模式，检出率，疾病类别 ，A类文献，B类文献,疾病描述] 增加[检出率（数字）、疾病类别（数字）]
    sort_Dic=defaultdict(list)
    sorted_Dic = {}
    for key,value in Dic.items():
        for line in value:
            ##获取检出率数字用于排序
            if '%' in line[3]:
                num_jianchu=re.search(ur'([\d\.]+)\\%',line[3]).group(1)
            elif '不明' in line[3]:
                num_jianchu='0'
            elif 'rare' in line[3]:
                num_jianchu='0.01'
            else:
                num_jianchu='101'
            line.append(num_jianchu)
            ##获取疾病种类转换为数字用于排序
            st=re.match('([AB])',line[4])
            if st.group(1)=='A':
                num_type='2'
            elif st.group(1)=='B':
                num_type='1'
            else:
                sys.stderr.write(u'请检查疾病种类')
                sys.exit(1)
            line.append(num_type)
            sort_Dic[key].append(line)
    for key in sort_Dic:
        sorted_Dic[key]=sorted(Dic[key],key=lambda x:(-float(x[8]),-float(x[9])))
    return sorted_Dic

def get_chr(Chr,start,end):
    site=''
    if start=='.' or start=='-' or Chr=='.' or Chr=='-':
        sys.stderr.write(u'染色体号或位点起始为空')
        sys.exit(1)
    if end=='.' or end=='-' or start == end:
        site='Chr'+Chr+':'+start
    else:
        site='Chr'+Chr+':'+start+' - '+end
    return site

def get_level_value(level):
    level_value=0
    if level=='致病': level_value=1
    elif level=='可能致病': level_value=2
    elif '临床意义未明1级' in re.sub(' ','',level): level_value=3
    elif '临床意义未明2级' in re.sub(' ','',level): level_value=4
    elif '临床意义未明3级' in re.sub(' ','',level): level_value=5
    elif '临床意义未明4级' in re.sub(' ','',level): level_value=6
    elif '临床意义未明5级' in re.sub(' ','',level): level_value=7
    elif level=='可能良性': level_value=8
    elif level=='良性': level_value=9
    else:
        sys.stderr.write(u'请检查致病等级')
        sys.exit(1)
    return level_value

def get_mut_type(mute):
    mute_Dic={'frameshift deletion':'移码缺失变异',
    'frameshift insertion':'移码插入变异',
    'frameshift substitution':'移码替换变异',
    'nonframeshift deletion':'非移码缺失变异',
    'nonframeshift insertion':'非移码插入变异',
    'nonframeshift substitution':'非移码替换变异',
    'stopgain':'无义变异',
    'stoploss':'终止缺失变异',
    'synonymous SNV':'同义变异',
    'nonsynonymous SNV':'错义变异',
    'splicing SNV':'剪切区单核苷酸变异',
    'splicing INDEL':'剪切区插入缺失变异',
    'intronic':'内含子变异',
    'SNV':'单核苷酸变异',
    }
    if mute_Dic.has_key(mute):
        mute_des=mute_Dic[mute]
    else:
        sys.stderr.write(u'请检查filter文件第9列是否为%s'%mute_Dic.keys())
        sys.exit(1)
    return mute_des

def read_CSdata():
    data='/PUBLIC/pipline/script/Report/XKFW/V3_180824/database/cs.data4.txt'
    disease_Dic=defaultdict(list)
    if os.path.isfile(data) == False:
        sys.stderr.write('%s is not a file' % data)
        sys.exit(1)
    try:
        with open(data, 'rU') as file_handle:
            for line in file_handle:
                if re.match('^疾病大类',line.strip()):continue
                line_1=re.sub('%','\%',line)
                #line_2=re.sub(r'\',r'\\',line_1)
                line_2=re.sub('_','\_',line_1)
                line_list=[x.strip() for x in line_2.strip().split('\t')]
                if len(line_list) < 9:
                    sys.stderr.write(u'cs.data.txt 列数不足 ')
		    sys.stderr.write(line_2)
                    sys.exit(1)
                disease_key=re.sub(' ','',line_list[1])+':'+re.sub(' ','',line_list[3])       #gene:disease_subtype作为key
                disease_Dic[disease_key]=[line_list[0],line_list[3],line_list[4],line_list[5],line_list[2],line_list[7],line_list[8]]  # 疾病名，疾病亚型，遗传模式，检出率，疾病类别 ，A类文献，B类文献
    except Exception as e:
        traceback.print_exc(e)
    return disease_Dic

def get_zongshu(site_Dic):
    cover_Dic={}
    hom_gene,hom_level,het_gene,het_level = tongji_hethom(site_Dic)
    hom_gene_set=[x for x in set(hom_gene)]
    het_gene_set=[x for x in set(het_gene)]
    if len(het_gene)==0:    ##全为纯合突变
        if len(hom_gene)==1:
            zongshu='本次受检样本中共检出 1 个纯合变异，存在于%s基因上，'%(hom_gene[0])
        elif len(set(hom_gene))==1:
            zongshu='本次受检样本中共检出 %s 个纯合变异，均存在于%s基因上，'%(str(len(hom_gene)),hom_gene_set[0])
        else:
            gene='、'.join(hom_gene)
            zongshu='本次受检样本中共检出 %s 个纯合变异，分别存在于%s基因上，'%(str(len(hom_gene)),gene)
    elif len(hom_gene)==0:  ##全为杂合突变
        if len(het_gene)==1:
            zongshu='本次受检样本中共检出 1 个杂合变异，存在于%s基因上，'%(het_gene[0])
        elif len(set(het_gene))==1:
            zongshu='本次受检样本中共检出 %s 个杂合变异，均存在于%s基因上，'%(str(len(het_gene)),het_gene_set[0])
        else:
            gene='、'.join(het_gene)
            zongshu='本次受检样本中共检出 %s 个杂合变异，分别存在于%s基因上，'%(str(len(het_gene)),gene)
    else:
        zongshu='本次受检样本中共检出 %s 个变异，%s个杂合变异和%s个纯合变异，分别存在于%s和%s基因上，'%(str(len(het_gene)+len(hom_gene)),str(len(het_gene)),str(len(hom_gene)),'、'.join(het_gene),'、'.join(hom_gene))
    zongshu += get_level(hom_level,het_level)
    zongshu += '与上述基因相关的疾病及检出位点详细信息请见下表。'
    cover_Dic['zongshu'] = zongshu
    return cover_Dic

def tongji_hethom(site_Dic):
    hom_gene = []
    het_gene = []
    hom_level = []
    het_level = []
    for key in sorted(site_Dic,key=lambda x:site_Dic[x][0][10]):
        value=site_Dic[key]
        for line in value:
            if line[4] == "纯合":
                hom_gene.append(key)
                hom_level.append(line[6])
            else:
                het_gene.append(key)
                het_level.append(line[6])
    return hom_gene,hom_level,het_gene,het_level

def get_level(hom,het):
    zong=hom+het
    Filt=set(zong)
    filt=[i for i in Filt]
    miaoshu=''
    if len(zong)==1:
        miaoshu = '严格遵循 ACMG 解读标准，判定为%s。'%zong[0]
    elif len(filt)==1:
        miaoshu = '严格遵循 ACMG 解读标准，均判定为%s。'%filt[0]
    elif len(het)==0:
        hom_miaoshu='、'.join(hom)
        miaoshu = '严格遵循 ACMG 解读标准，依次判定为%s。'%hom_miaoshu
    elif len(hom)==0:
        het_miaoshu='、'.join(het)
        miaoshu = '严格遵循 ACMG 解读标准，依次判定为%s。'%het_miaoshu
    else:
        hom_miaoshu='、'.join(hom)
        het_miaoshu='、'.join(het)
        miaoshu = '严格遵循 ACMG 解读标准，依次判定为%s和%s。'%(het_miaoshu,hom_miaoshu)
    return miaoshu

def get_disease_table(Dic,siteDic,BZDic):
    #疾病名0，疾病亚型1，遗传模式2，检出率3，疾病类别 4，A类文献5，B类文献6,疾病描述7，检出率数字8，疾病类别数字9
    cover_Dic={}
    table=''
    mode,litea,liteb = prepare_zhushi(Dic)
    Dic_new=Dic
    sorted_Dic=siteDic
    miaoshubz=''
    flag=0;
    for key in sorted(sorted_Dic,key=lambda x:sorted_Dic[x][0][10]):  ##根据致病性排序
        if len(Dic_new[key]) ==1:
            table += r'''
%s & %s & %s & %s & %s \\ \hline'''%(key,Dic_new[key][0][1],Dic_new[key][0][2],Dic_new[key][0][3],Dic_new[key][0][4])
        else:
            table += get_rows(key,Dic_new[key])
        if key in BZDic:
            if flag==0:
                  miaoshubz='注：'
            flag=1 
            joina='、'.join(BZDic[key])
            miaoshubz+='%s基因可能还与%s相关，'%(key,joina)
    if flag==1:
         miaoshubz+='但目前研究较少，如已经出现该病的临床表型，请结合临床考虑该位点的致病性。\\'
    cover_Dic['result_table']=table
    cover_Dic['result_zhujie']=miaoshubz
    cover_Dic['genetic_model']=get_model(mode)
    cover_Dic['typea_num'],cover_Dic['typeb_num'],cover_Dic['type_lite']=get_lite(litea,liteb)
    return cover_Dic

def get_lite(lia,lib):
    lite=''
    liteanum=''
    litebnum=''
    lia_New = prepare_lite(lia)
    lib_New = prepare_lite(lib)
    lia_new = [i for i in set(lia_New)]
    lib_new = [i for i in set(lib_New)]
    li=lia_new + lib_new
    ##获取文献列表
    if len(li)==0:
        pass
    else:
        lite=r'''
\par \zihao{6}{\color{DarkBlue} 参考文献}\\ \vspace*{-6mm}
\begin{spacing}{1.5}
\zihao{-6} \color{MyFontGray}
\setdefaultleftmargin{2em}{}{}{}{.5em}{.5em}
\begin{compactenum}'''
        for index,item in enumerate(li):
            lite += r'''
\item %s'''%item
        lite += r'''
\end{compactenum}
\end{spacing}'''
    ##获取索引数字
    suoyin=[]
    if len(lia_new)==0:
        liteanum=''
    else:
        for index,item in enumerate(lia_new):
            suoyin.append(str(index+1))
        liteanum = ','.join(suoyin)
        liteanum = '[' + liteanum +']'
    suoyin=[]
    if len(lib_new)==0:
        litebnum=''
    else:
        for index,item in enumerate(lib_new):
            suoyin.append(str(index+1+len(lia_new)))
        litebnum = ','.join(suoyin)
        litebnum = '[' + litebnum +']'
    return liteanum,litebnum,lite

def prepare_lite(li):
    li_new=[]
    for i in li:
        if '|' in i:
            arr=i.split('|')
            for jj in arr:
	        if jj=='' or jj =='-' or jj =='.':
                    pass
                else:
                    li_new.append(jj)	
                #li_new += arr
        else:
            li_new.append(i)
    return li_new

def get_model(li):
    genetic_model=''
    mode_Dic={'AD':'AD：常染色体显性遗传；',
    'AR':'AR：常染色体隐性遗传；',
    'XLD':'XLD：X 染色体连锁显性遗传；',
    'XLR':'XLR：X 染色体连锁隐性遗传；',
    'XL':'XL：X 染色体连锁遗传；',
    'AD/AR':'AD/AR：常染色体遗传；',
    'AR/AD':'AR/AD：常染色体遗传；'}
    for mode in li:
        if mode_Dic.has_key(mode):
            genetic_model += mode_Dic[mode]
        elif mode == r'不明':
            pass
        else:
            sys.stderr.write(u'遗传模式为%s，不在AD AR XLD XLR XL AD/AR AR/AD内'%mode)
            sys.exit(1)
    return genetic_model

def get_rows(gene,li):
    m=0
    table=''
    for line in li:
        m += 1
        if m==1:
            table += r'''
\multirow{%s}{*}{%s} & %s & %s & %s & %s \\ \cline{2-5}'''%(str(len(li)),gene,line[1],line[2],line[3],line[4])
        elif m==len(li):
            table += r'''
& %s & %s & %s & %s \\ \hline'''%(line[1],line[2],line[3],line[4])
        else:
            table +=r'''
& %s & %s & %s & %s \\ \cline{2-5}'''%(line[1],line[2],line[3],line[4])
    return table

def prepare_zhushi(Dic):
    # [疾病名，疾病亚型，遗传模式，检出率，疾病类别 ，A类文献，B类文献,疾病描述,检出率（数字）、疾病类别（数字）
    literatureA= []
    literatureB= []
    genetic_mode=[]
    for key,value in Dic.items():
        for line in value:
            genetic_mode.append(line[2])   ##获取遗传模式
            if line[5] !='-' and line[5] !='' :literatureA.append(line[5])   ##获取A类文献
            if line[6] !='-' and line[6] !='' :literatureB.append(line[6])   ##获取B类文献
    mode_New = set(genetic_mode)
    litea_New = set(literatureA)
    liteb_New = set(literatureB)
    mode_new = [i for i in mode_New]
    litea_new = [i for i in litea_New]
    liteb_new = [i for i in liteb_New]
    return mode_new,litea_new,liteb_new

def get_site_table(Dic):
    table=''
    cover_Dic={}
    num=0
    sortDic=Dic
    #nue,acid,chr_site,exon,het/hom/,mute_type,zhibingxing,literature,gene_des,site_des
    for key in sorted(sortDic,key=lambda x:sortDic[x][0][10]):
        if len(sortDic[key])==1:
            color=get_color(sortDic[key][0][6])
            table += r'''
\color{MyFontGray}
%s  & %s & %s & %s &   \makecell*[c]{%s} & %s & %s & \textcolor{%s}{%s} \\ \hline
'''%(key,sortDic[key][0][0],sortDic[key][0][1],sortDic[key][0][2],sortDic[key][0][3],sortDic[key][0][4],sortDic[key][0][5],color,sortDic[key][0][6],)
        else:
            num=get_multline(sortDic[key])
            m=0
            for line in sortDic[key]:
                m += 1
                color=get_color(line[6])
                if m==1:
                    table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s} & %s & %s & %s &   \makecell*[c]{%s} & %s & %s & \textcolor{%s}{%s} \\ \cline{2-8}'''%(str(num),key,line[0],line[1],line[2],line[3],line[4],line[5],color,line[6])
                elif m==len(Dic[key]):
                    table += r'''
\color{MyFontGray}
& %s & %s & %s &   \makecell*[c]{%s} &  %s & %s & \textcolor{%s}{%s} \\ \hline'''%(line[0],line[1],line[2],line[3],line[4],line[5],color,line[6])
                else:
                    table += r'''
\color{MyFontGray}
& %s & %s & %s &   \makecell*[c]{%s} &  %s & %s & \textcolor{%s}{%s} \\ \cline{2-8}'''%(line[0],line[1],line[2],line[3],line[4],line[5],color,line[6])
    cover_Dic['site_table']=table
    return cover_Dic

def get_fubiao_table(Dic,siteDic):
    cover_Dic={}
    table=''
    num=0
    sortDic=siteDic
    wenxian=''
    DicK=[]
    for key in sorted(sortDic,key=lambda x:sortDic[x][0][8]):
        m=0
        m1=len(sortDic[key])
        m2=len(Dic[key])
        print key,m1,m2
        if len(Dic[key])==1:
            for line in sortDic[key]:
                table += r'''
\color{MyFontGray}
%s & %s & %s &   \makecell*[c]{%s} & %s & %s & %s & \makecell*[c]{%s} & %s & %s \\ \hline'''%(key,line[0],line[1],line[2],line[3],line[4],line[6],line[7],Dic[key][0][1],Dic[key][0][2])
        elif len(Dic[key])>1:
            num=len(Dic[key])
            for line in sortDic[key]:
                mm=0
                for line2 in Dic[key]:
                    mm+=1
                    if mm==1:
                        strints=line[7]
                        strintsnew=re.sub('临床意义未','临床意义\\\\\\\\未',strints)
                        value=len(line[4])
                        print value
                        print line[4]
                        if value <= 12 :
                            table += r'''
\color{MyFontGray}
\raisebox{-2mm}{\multirow{%s}{*}{%s}}  & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} &  %s & %s \\ \cline{9-10}'''%(str(num),key,str(num),line[0],str(num),line[1],str(num),line[2],str(num),line[3],str(num),line[4],str(num),line[6],str(num),strintsnew,line2[1],line2[2])
                        else:
                            bianyitype=line[4][0:12]+'\\\\'+line[4][12:value]
                            print bianyitype
                            table += r'''
\color{MyFontGray}
\raisebox{-2mm}{\multirow{%s}{*}{%s}}  & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} &  %s & %s \\ \cline{9-10}'''%(str(num),key,str(num),line[0],str(num),line[1],str(num),line[2],str(num),line[3],str(num),bianyitype,str(num),line[6],str(num),strintsnew,line2[1],line2[2])
                    elif mm==len(Dic[key]):
                        table += r'''
\color{MyFontGray}
&   &  &    &   &  &  &  & %s & %s  \\ \hline'''%(line2[1],line2[2])
                    else:
                        table += r'''
\color{MyFontGray}
&   &  &    &   &  & &  & %s & %s   \\ \cline{9-10}'''%(line2[1],line2[2])
        DicK.append(Dic[key])
    wenxian = get_wenxianfu(DicK)
    cover_Dic['fubiao_table']=table
    cover_Dic['fubiao_wenxian']=wenxian
    return cover_Dic

def get_fubiao_table_bak(Dic,siteDic):
    cover_Dic={}
    table=''
    num=0
    sortDic=siteDic
    wenxian=''
    DicK=[]
    for key in sorted(sortDic,key=lambda x:sortDic[x][0][8]):
    #for key in  sortDic:
        m=0
        m1=len(sortDic[key])
        m2=len(Dic[key])
        print key,m1,m2
        if len(sortDic[key])==1 and len(Dic[key])==1:
            table += r'''
\color{MyFontGray}
%s  & %s & %s &  \makecell*[c]{%s} & %s & %s & %s & %s & %s  & %s \\ \hline
'''%(key,sortDic[key][0][0],sortDic[key][0][1],sortDic[key][0][2],sortDic[key][0][3],sortDic[key][0][4],sortDic[key][0][6],sortDic[key][0][7],Dic[key][0][1],Dic[key][0][2])
        elif len(sortDic[key]) > 1 and len(Dic[key])==1:
            num=len(sortDic[key])
            for line in sortDic[key]:
                m += 1
                if m==1:
                    table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s}  & %s & %s &   \makecell*[c]{%s} & %s & %s & %s & %s & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} \\ \cline{2-8}'''%(str(num),key,line[0],line[1],line[2],line[3],line[4],line[6],line[7],str(num),Dic[key][0][1],str(num),Dic[key][0][2])
                elif m==num:
                    table += r'''
\color{MyFontGray}
& %s  & %s &   \makecell*[c]{%s} &  %s & %s & %s  & %s  & &\\ \hline'''%(line[0],line[1],line[2],line[3],line[4],line[6],line[7])
                else:
                    table += r'''
\color{MyFontGray}
& %s & %s &   \makecell*[c]{%s} &  %s & %s & %s  & %s & & \\ \cline{2-8}'''%(line[0],line[1],line[2],line[3],line[4],line[6],line[7])#wph 180709
        elif m1 == 1 and  m2 > 1:
            mm=0
            num=len(Dic[key])
            for line in Dic[key]:
                mm+=1
                if mm==1:
                    strints=sortDic[key][0][7]
                    strintsnew=re.sub('临床意义未','临床意义\\\\\\\\未',strints)
                    value=len(sortDic[key][0][4])
                    print value
                    print sortDic[key][0][4]
                    if value <= 12 :
                         table += r'''
\color{MyFontGray}
\raisebox{-2mm}{\multirow{%s}{*}{%s}}  & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} &  %s & %s \\ \cline{9-10}'''%(str(num),key,str(num),sortDic[key][0][0],str(num),sortDic[key][0][1],str(num),sortDic[key][0][2],str(num),sortDic[key][0][3],str(num),sortDic[key][0][4],str(num),sortDic[key][0][6],str(num),strintsnew,line[1],line[2])           
                    else:
                        bianyitype=sortDic[key][0][4][0:12]+'\\\\'+sortDic[key][0][4][12:value]
                        print bianyitype
                        table += r'''
\color{MyFontGray}
\raisebox{-2mm}{\multirow{%s}{*}{%s}}  & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} &  %s & %s \\ \cline{9-10}'''%(str(num),key,str(num),sortDic[key][0][0],str(num),sortDic[key][0][1],str(num),sortDic[key][0][2],str(num),sortDic[key][0][3],str(num),bianyitype,str(num),sortDic[key][0][6],str(num),strintsnew,line[1],line[2])  
                elif mm==len(Dic[key]):
                    table += r'''
\color{MyFontGray}
&   &  &    &   &  &  &  & %s & %s  \\ \hline'''%(line[1],line[2]) 
                else:
                    table += r'''
\color{MyFontGray}
&   &  &    &   &  & &  & %s & %s   \\ \cline{9-10}'''%(line[1],line[2])      
        else:
            if m1 == m2:
                num=len(Dic[key])
                for index in range(0,num):
                    line2=Dic[key][index]
                    line1=sortDic[key][index]
                    if index==0:
                        table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s} &  %s & %s  &  \makecell*[c]{%s}  &  %s & %s  & %s  & %s & %s & %s  \\ \cline{2-10}'''%(str(num),key,line1[0],line1[1],line1[2],line1[3],line1[4],line1[6],line1[7],line2[1],line2[2]) 
                    elif index==num-1: 
                        table += r'''
&  %s & %s  &  \makecell*[c]{%s}  &  %s & %s  & %s  & %s & %s & %s \\ \hline'''%(line1[0],line1[1],line1[2],line1[3],line1[4],line1[6],line1[7],line2[1],line2[2]) 
		    else:
		        table += r'''
&  %s & %s  &  \makecell*[c]{%s}  &  %s & %s  & %s  & %s & %s & %s \\ \cline{2-10}'''%(line1[0],line1[1],line1[2],line1[3],line1[4],line1[6],line1[7],line2[1],line2[2])
            elif m2==2: 
                pass
	    elif m1>m2 and m2>1:
		diseasename=""
		hercent=""
		for line in Dic[key]:
		    diseasename=diseasename+"\\\\"+line[1]
                    hercent=hercent+"\\\\"+line[2]
                diseasename=diseasename.strip("\\")
                hercent=hercent.strip("\\")
                m=0;
                num=len(sortDic[key])
                for line in sortDic[key]:
                    m += 1
                    if m==1:
                        table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s}  & %s & %s &   \makecell*[c]{%s} & %s & %s & %s & %s & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} \\ \cline{2-8}'''%(str(num),key,line[0],line[1],line[2],line[3],line[4],line[6],line[7],str(num),diseasename,str(num),hercent)
                    elif m==num:
                        table += r'''
\color{MyFontGray}
& %s  & %s &   \makecell*[c]{%s} &  %s & %s & %s  & %s  & &\\ \hline'''%(line[0],line[1],line[2],line[3],line[4],line[6],line[7])
                    else:
                        table += r'''
\color{MyFontGray}
& %s & %s &   \makecell*[c]{%s} &  %s & %s & %s  & %s & & \\ \cline{2-8}'''%(line[0],line[1],line[2],line[3],line[4],line[6],line[7])
            else:
                m11=m2/2
                m22=m2-m11
                ac=""
                aa=""
                tr=""
                av=""
                avc=""
                achr=""
                achrg=""#wph add 180703
                mm=0
                for line in sortDic[key]:
                    mm+=1
                    if mm>1:
                        ac=ac+"\\\\"+line[0]
                        aa=aa+"\\\\"+line[1]
                        tr=tr+"\\\\"+line[2]
                        av=av+"\\\\"+line[3]
                        avc=avc+"\\\\"+line[4]
                        achr=achr+"\\\\"+line[6]
                        achrg=line[7]
                ac=ac.strip("\\")
                aa=aa.strip("\\")
                tr=tr.strip("\\")
                av=av.strip("\\")
                avc=avc.strip("\\")
                achr=achr.strip("\\")
                achrg=achrg.strip("\\")
                mm=0
                num=len(Dic[key])
                for line in Dic[key]:
                    mm+=1
                    if mm==1 and m11==1:
                        strints=sortDic[key][0][7]
                        strintsnew=re.sub('临床意义未','临床意义\\\\\\\\未',strints)
                        table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s}  & %s & %s & \makecell*[c]{%s} & %s & %s & %s & \makecell*[c]{%s}  & %s  & %s  \\ \cline{2-10}'''%(str(num),key,sortDic[key][0][0],sortDic[key][0][1],sortDic[key][0][2],sortDic[key][0][3],sortDic[key][0][4],sortDic[key][0][6],strintsnew,line[1],line[2])
                    elif mm==1 and m11>1:
                        strints=sortDic[key][0][7]
                        print strints
                        strintsnew=re.sub('临床意义未','临床意义\\\\\\\\未',strints)
                        value=len(sortDic[key][0][4])
                        print value
                        print sortDic[key][0][4]
                        if value <= 12 :
                            table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s}  & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & %s  & %s \\ \cline{8-10}'''%(str(num),key,m11,sortDic[key][0][0],m11,sortDic[key][0][1],m11,sortDic[key][0][2],m11,sortDic[key][0][3],m11,sortDic[key][0][4],m11,sortDic[key][0][6],m11,strintsnew,line[1],line[2])
                        else:
                            bianyitype=sortDic[key][0][4][0:12]+'\\\\'+sortDic[key][0][4][12:value]
                            print bianyitype
                            table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s}  & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{%s}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & %s  & %s \\ \cline{8-10}'''%(str(num),key,m11,sortDic[key][0][0],m11,sortDic[key][0][1],m11,sortDic[key][0][2],m11,sortDic[key][0][3],m11,bianyitype,m11,sortDic[key][0][6],m11,strintsnew,line[1],line[2])
                    elif mm ==m11:
                        table += r'''
\color{MyFontGray}
  &  &  &  &  &  &  &  & %s  & %s  \\ \cline{2-8}\cline{9-10}'''%(line[1],line[2])
                    elif mm== m11+1 and m22==1:
                        table += r'''
\color{MyFontGray}
  & %s & %s & \makecell*[c]{%s} & %s & %s & %s & %s  & %s  & %s  \\ \cline{9-10}'''%(ac,aa,tr,av,avc,achr,sortDic[key][0][7],line[1],line[2])   
                    elif mm== m11+1 and m22>1:
                        strints=achrg
                        strintsnew=re.sub('临床意义未','临床意义\\\\\\\\未',strints)
                        table += r'''
\color{MyFontGray}
  & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}} & \raisebox{-2mm}{\multirow{%s}{*}{\makecell*[c]{%s}}}  & %s  & %s  \\ \cline{9-10}'''%(m22,ac,m22,aa,m22,tr,m22,av,m22,avc,m22,achr,m22,strintsnew,line[1],line[2])
                    elif mm==len(Dic[key]):
                        table += r'''
\color{MyFontGray}
&   &  &    &   &  &  &  & %s & %s  \\ \hline'''%(line[1],line[2])
                    else:
                        table += r'''
\color{MyFontGray}
&   &  &    &   &  & &  & %s & %s & %s  \\ \cline{9-10}'''%(line[1],line[2])    
        DicK.append(Dic[key])
    wenxian = get_wenxianfu(DicK)
    cover_Dic['fubiao_table']=table
    cover_Dic['fubiao_wenxian']=wenxian
    return cover_Dic


def get_duanzi(name):
    new_name=''
    if 'ins' in name:
        new_name=re.sub('ins',' ins',name)
    elif 'del' in name:
        new_name=re.sub('del',' del',name)
    else:
        new_name=name
    return new_name

def get_multline(li):
    level=[]
    num=0
    for line in li:
        if u'临床意义未明' in line[6]:
            num += 2
        else:
            num += 1
    return num

def get_color(cli):
    color=''
    if cli ==u"致病": color="red"
    elif cli ==u"可能致病": color="orange"
    elif cli ==u"良性" or cli ==u"可能良性": color="MyGreen"
    else: color="yellow"
    return color

def get_jiexi_table(disease,site):
    site_sortDic=site
    disease_Dic=disease
    cover_Dic={}
    jiexi=''
    table_li=[]
    jiexi_beizhu=r'''
\par \vspace*{2mm}
\color{DarkBlue}{\zihao{6} 备注}\\
\zihao{-6} \color{MyFontGray} 预测软件包括 SIFT、PolyPhen2、M-CAP 和 REVEL, 用于预测错义变异是否会导致蛋白结构和功能发生改变, 准确率大致在 65\%-80\%。需要说明的是, 各软件提供的只是基于算法的预测结果, 可以为位点解读提供参考, 但不足以作为致病性判定的唯一标准, 还需要结合其他信息综合判断。\\ \\
以上位点变异为基因检测的客观结果，但变异是否真正致病需要由医生结合临床信息综合判断。如果您想进一步了解基因变异的信息或家族遗传史 ，建议联系家人一起进行检测。'''
    jiexi_des_head=r'''
\vspace*{-4mm}
\begin{spacing}{1.5}
\color{MyFontGray} \zihao{-5}
\rowcolors{1}{}{LightBlue}
\arrayrulecolor{DarkBlue}
\renewcommand\arraystretch{1.5}
\begin{longtable}{C{3cm} L{13.2cm}}'''
    jiexi_des_foot=r'''
\hline
\end{longtable}
\end{spacing}'''
    for key in sorted(site_sortDic,key=lambda x:site_sortDic[x][0][10]):
        jiexi=''
        disease_li,disease_des,site_li,site_des,level=prepare_jiexi(site_sortDic[key],disease_Dic[key])
        if len(site_li)==1:
            jiexi +=r'''
\begin{spacing}{1.2}
\zihao{5} \color{white} \sym \sye
\renewcommand\arraystretch{1}
\colorbox{DarkBlue}{
\begin{tabular}{L{2cm} L{14cm}}
风险位点: & %s—— %s\\
相关疾病: & %s\\
位点致病性: & %s \\
\end{tabular}}
\end{spacing} '''%(key,site_li[0],'、'.join(disease_li),level[0])
            jiexi += jiexi_des_head
            jiexi += r'''
变异解析 & %s \\'''%site_des[0]
            jiexi += r'''
基因描述 & %s \\'''%site[key][0][8]
            jiexi += get_disease_row(disease_li,disease_des)
            jiexi += jiexi_des_foot
            jiexi += get_wenxian(site_sortDic[key])
            jiexi += jiexi_beizhu
            table_li.append(jiexi)
        else:
            jiexi += get_head(key,site_li,disease_li,level)
            jiexi += jiexi_des_head
            for index,item in enumerate(site_des):
                jiexi +=r'''
变异解析%s & %s \\'''%(str(index+1),item)
            jiexi += r'''
基因描述 & %s \\'''%site[key][0][8]
            jiexi += get_disease_row(disease_li,disease_des)
            jiexi += jiexi_des_foot
            jiexi += get_wenxian(site_sortDic[key])
            jiexi += jiexi_beizhu
            table_li.append(jiexi)
    table=r'\newpage '.join(table_li)
    cover_Dic['jiexi_table']=table
    return cover_Dic

def get_wenxian(li):
    wenxian=[]
    lite=''
    for i in li:
        if i[7]=='' or i[7] =='-':
            pass
        else:
            wenxian.append(i[7])
    wenxian_new=prepare_lite(wenxian)
    if len(wenxian_new) ==0:
        lite=r'''
\par \vspace*{5mm}'''
    else:
        lite =r'''
\par \vspace*{5mm}
\zihao{6}{\color{DarkBlue} 参考文献}\\ \vspace*{-6mm}
\begin{spacing}{1.5}
\zihao{-6} \color{MyFontGray}
\setdefaultleftmargin{2em}{}{}{}{.5em}{.5em}
\begin{compactenum}'''
        for i in set(wenxian_new):
            lite +=r'''
\item %s'''%i
        lite += r'''
\end{compactenum}
\end{spacing}'''
    return lite


def get_wenxianfu(li):
    wenxian=[]
    lite=''
    for i in li:
        for j in i:
            if j[5]=='' or j[5] =='-' or j[5] =='.':
                pass
            else:
                wenxian.append(j[5])
            if j[6]=='' or j[6] =='-'  or j[6] =='.':
                pass
            else:
                wenxian.append(j[6])
    wenxian_new=prepare_lite(wenxian)
    if len(wenxian_new) ==0:
        lite=r'''
\par \vspace*{5mm}'''
    else:
        lite =r'''
\par \vspace*{5mm}
\zihao{6}{\color{DarkBlue} 参考文献}\\ \vspace*{-6mm}
\begin{spacing}{1.5}
\zihao{-6} \color{MyFontGray}
\setdefaultleftmargin{2em}{}{}{}{.5em}{.5em}
\begin{compactenum}'''
        for i in set(wenxian_new):
            lite +=r'''
\item %s'''%i
        lite += r'''
\end{compactenum}
\end{spacing}'''
    return lite


def get_head(gene,site,disease,level):
    head=''
    if len(site) != len(level):
        sys.stderr.write(u'位点数与致病等级数不一致')
        sys.exit(1)
    else:
        head =r'''
\begin{spacing}{1.2}
\zihao{5} \color{white} \sym \sye
\renewcommand\arraystretch{1}
\colorbox{DarkBlue}{
\begin{tabular}{L{2cm} L{14cm}}'''
        for index,value in enumerate(site):
            head +=r'''
风险位点%s: & %s—— %s, %s\\'''%(str(index+1),gene,value,level[index])
        head +=r'''
相关疾病: & %s\\
\end{tabular}}
\end{spacing} '''%('、'.join(disease))
    return head

def get_disease_row(disease,disease_des):
    table=''
    if len(disease) != len(disease_des):
        sys.stderr.write(u'请检查filter是否每个位点都包含疾病名称、疾病描述')
        sys.exit(1)
    else:
        for index,item in enumerate(disease):
            table += r'''
%s & %s \\'''%(item,disease_des[index])
    return table


def prepare_jiexi(site_li,disease_li):
    #nue,acid,chr_site,exon,het/hom/,mute_type,zhibingxing,literature,gene_des,site_des
    disease=[]
    disease_des=[]
    site=[]
    site_des=[]
    level=[]
    for line in disease_li:
        if line[0] not in disease:
            disease.append(line[0])
            disease_des.append(line[7])
        else:
            pass
    for line in site_li:
        site_detail=line[0]+':'+line[1]+', '+line[4]
        site.append(site_detail)
        site_des.append(line[9])
        level.append(line[6])
    return disease,disease_des,site,site_des,level


def main():
    argInit()
    cover_Dic={}
    cover_Dic['part'] = part
    if module.upper() == 'CS':
        cover_Dic['title']='Part'+part+' — 检测结果'
    else:
        cover_Dic['title']='Part'+part+' — 遗传性心血管病基因检测结果'
    out_CS=open(outDir+'/'+sampleid+'_Cardio_CS.tex','w')
    if os.path.isfile(yinyang_file) == False:
        sys.stderr.write('%s is not a file' % yinyang_file)
        sys.exit(1)
    try:
        yinyang_handle = open(yinyang_file,'r')
        yinyang = yinyang_handle.readlines()[0].strip()
        yinyang_handle.close()
        if yinyang == '0':
            id_filter,id_filterfu,id_filterbz = get_filter(filterFile,sampleid)
            with open('/PUBLIC/pipline/script/Report/XKFW/V3_180824/tempt/Cardio_CSyin_temple_2.tex','r') as file_handle:
                CS_temp = Template(file_handle.read())
                diseaseDic,siteDic,diseaseDic_fu,siteDic_fu,miaoshubzz  = prepare_allDic(id_filter,id_filterfu,id_filterbz,0)
                fubiao_table = get_fubiao_table(diseaseDic_fu,siteDic_fu)
                cover_Dic.update(fubiao_table)
                num=len(id_filterfu)
                if num==0:
                    cover_Dic['yin']=r'''此次遗传性心血管病基因检测，包含179个疾病相关基因，筛查36种疾病亚型。检测结果显示您并未发生明确致病变异的位点，因此，您的疾病症状并不是由这些已知基因的致病变异引起的，可能是目前研究尚未涉及的基因变异，也可能是由于其他因素，如不良生活习惯、环境、其他疾病等引起。建议您注意生活细节，积极改变不良生活方式，并配合医生治疗。 

以上检测结果仅基于当前的科学研究，随研究文献的更新，疾病-基因的对应关系可能会发生变化，因此我们的报告具有一定的时效性 。
'''
                else:
                    cover_Dic['yin']=r'''此次遗传性心血管病基因检测，包含179个疾病相关基因，筛查36种疾病亚型。检测结果显示您并未发生高风险的位点变异；一些变异位
点由于相关研究较少、ACMG证据条目不足，基于现有研究判断其致病性等级较低，但不排除导致疾病的可能性，这类位点的致病性等级及相关证据条目展示在附表中，供您和临床医生参考。

若检出位点的相关疾病与您现有的临床表型不符，那么您的相关临床症状可能与其他未知基因变异、疾病的不完全外显性、不良生活习惯或环境等因素相关。建议您注意生活细节，积极改变
不良生活方式，或去医院进行相关疾病的进一步检测，听从医生的临床建议，定期体检，并采取积极的生活干预。

若您未发现任何相关疾病表型，则这类 ACMG 致病类证据较少的位点可能并不致病。 

以上检测结果仅基于当前的科学研究，随研究文献的更新，位点致病性情况可能会发生变化，属于正常现象。如果您想进一步了解基因变异的信息或家族遗传史，建议联系家人一起进行检测
。
'''

                report = CS_temp.safe_substitute(cover_Dic)
                print >> out_CS,report
        else:
            id_filter,id_filterfu,id_filterbz = get_filter(filterFile,sampleid)
            with open('/PUBLIC/pipline/script/Report/XKFW/V3_180824/tempt/Cardio_CSyang_temple.tex','r') as file_handle:
                CS_temp = Template(file_handle.read())
                diseaseDic,siteDic,diseaseDic_fu,siteDic_fu,miaoshubzz = prepare_allDic(id_filter,id_filterfu,id_filterbz,1)
                zongshu = get_zongshu(siteDic)
                disease_table = get_disease_table(diseaseDic,siteDic,miaoshubzz)
                site_table = get_site_table(siteDic)
                fubiao_table = get_fubiao_table(diseaseDic_fu,siteDic_fu)
                jiexi_table = get_jiexi_table(diseaseDic,siteDic)
                cover_Dic.update(zongshu)
                cover_Dic.update(disease_table)
                cover_Dic.update(site_table)
                cover_Dic.update(fubiao_table)
                cover_Dic.update(jiexi_table)
                report = CS_temp.safe_substitute(cover_Dic)
                print >> out_CS,report
    except Exception as e:
        traceback.print_exc(e)
    finally:
        out_CS.close()

if __name__ == '__main__':
    main()
