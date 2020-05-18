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
        '-site', help='file GXY.filter.txt', required=True)
    parser.add_argument(
        '-part', help='2 / 3 / 4', required=True)
    parser.add_argument(
        '-out', help='project name', default=os.getcwd())
    parser.add_argument(
        '-sampleID', help='sampleid', required=True)
    parser.add_argument(
        '-disease', help='GXY.disease.table', required=True)
    parser.add_argument(
        '-module', help='module', required=True)
    parser.add_argument(
        '-yinyang', help='阴阳性', required=True)
    argv = vars(parser.parse_args())
    global filterFile, part,outDir,sampleid,disease_file,module,yinyang_file
    filterFile = argv['site'].strip()
    part = argv['part'].strip()
    outDir = os.path.abspath(argv['out'].strip())
    sampleid=argv['sampleID'].strip()
    disease_file=argv['disease'].strip()
    module=argv['module'].strip()
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
                if re.match('^Sample',line.strip()):continue
                if re.match('^####',line.strip()):continue
                line_list=[x.strip() for x in line.strip().split('\t')]
                if re.match('^###',line.strip()):
                    newsam=re.sub('###','',line_list[0])
                    if  newsam==sampleid:
                        line_1=re.sub('%','\%',line)
                        line_2=re.sub('_','\_',line_1)
			line_3=re.sub('~','-',line_2) ### panqi 20181213
                        filter_bz.append(line_3.strip())
		elif re.match('^##N',line.strip()):
                    newsam=re.sub('##','',line_list[0])
                    if  newsam==sampleid:
                        line_1=re.sub('%','\%',line)
                        line_2=re.sub('_','\_',line_1)
			line_3=re.sub('~','-',line_2) ### panqi 20181213
                        filter_fu.append(line_3.strip())
		else:
                    id_li.append(line_list[0])
                    if ''.join(line_list[1:6]) == '-----':continue
                    if line_list[0]==sampleid:
                        line_1=re.sub('%','\%',line)
                        #line_2=re.sub('\\',r'$\backslash$',line_1)
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
    return filter_li,filter_fu

def prepare_allDic(files,disease_file,filesfu):
    disease_Dic = defaultdict(list)
    site_Dic = defaultdict(list)
    site_li = []
    gene_li= []
    key_li=[]
    het_hom ={'het':u'杂合','hom':u'纯合','unknown':u'未明','hem':u'半合子'} #wph add hem 181212
    data_Dic=read_GXYdata(disease_file)
    try:
        for line in files:
            line_list=[x.strip() for x in line.strip().split('\t')]
            if len(line_list) < 89:
                sys.stderr.write('filter文件列数不足，请检查')
                sys.exit(1)
            ##获取disease_Dic
            gene=line_list[85].strip()
            subtype_disease=re.sub(' ','',line_list[83].strip())
            key=gene+':'+subtype_disease
            if  key not in key_li:
                key_li.append(key)
                if data_Dic.has_key(key):
                    disease_detail=data_Dic[key]
                    biaoxian=line_list[89]
                    genetic_mode=line_list[84]
                    dalei_miaoshu=line_list[79]
                    zhibingxing=line_list[78]
                    zhiliao=line_list[88]
                    disease_detail.append(genetic_mode)
                    disease_detail.append(dalei_miaoshu)
                    disease_detail.append(biaoxian)
                    disease_detail.append(zhibingxing)
                    disease_detail.append(zhiliao)
                    disease_Dic[gene].append(disease_detail)  #key>gene：subtype_disease,value>[大类+亚类，亚类，大类，遗传模式，疾病描述，临床表现,致病性,治疗建议]
                else:
                    sys.stderr.write(u'请检查filter内%s与GXY.disease.table中的gene名、疾病亚型名是否一致'%key)#gene_dis_key2)
                    sys.exit(1)
            else:
                pass
            ##获取site_Dic
            site=line_list[1]+":"+line_list[2]+":"+line_list[3]+":"+line_list[4]+":"+line_list[5] #chr:start:end:ref:alt
            if site not in site_li:
                if line_list[90]=='-' or line_list[90]=='.' or line_list[90]=='.//.,':    ### panqi 180817
                    exon='-'
                else:
		    #ncbi_tran=line_list[12].split('//')   #changed by sunxq
                    exons=line_list[90].split(':')
		    #if len(ncbi_tran)>1 and ncbi_tran[1] != '' and  ncbi_tran[1] != '.,' :
		    #   exons=ncbi_tran[1].split(':')
		    #else:
                    #   exons=line_list[12].split(':')
                    exon=exons[0]+'\\\\'+exons[1]
                chr_site=get_chr(line_list[1],line_list[2],line_list[3])
                mute_type=get_mut_type(line_list[9])
                level_value=get_level_value(line_list[78])
                site_detail=[line_list[86],line_list[87],chr_site,exon,het_hom[line_list[31].lower()],mute_type,line_list[78],line_list[82],line_list[80],line_list[81],line_list[88],level_value,line_list[44],line_list[83],line_list[84]]
                                      #nue,     acid,    chr_site,exon,het/hom,                       mute_type,zhibingxing,    literature,  gene_des,     site_des,      advise,      level_value   true_false    subtype       genetic_model
                site_Dic[line_list[85]].append(site_detail)
                site_li.append(site)
            else:
                pass
    except Exception as e:
            traceback.print_exc(e)
    site_Dic_fubiao = defaultdict(list)
    site_li_fubiao = []
    try:
        for line in filesfu:
            line_list=[x.strip() for x in line.strip().split('\t')]
            site=line_list[1]+":"+line_list[2]+":"+line_list[3]+":"+line_list[4]+":"+line_list[5]+line_list[83]
            if site not in site_li_fubiao:
                if line_list[90]=='-' or line_list[90]=='.' or line_list[90]=='.//.,':
                    exon='-'
                else:
                    exons=line_list[90].split(':')
                    exon=exons[0]+'\\\\'+exons[1]
                    #ncbi_tran=line_list[12].split('//')
                    #if ncbi_tran[1] != ''  and ncbi_tran[1] != '.' and ncbi_tran[1] != '.,':
                    #    exons=ncbi_tran[1].split(':')
                    #else:
                    #    exons=line_list[12].split(':')
                    #exon=exons[1]+'\\\\'+exons[2]
                mute_type=get_mut_type(line_list[9])
                level_value=get_level_value(line_list[78])
                site_detail=[line_list[86],line_list[87],exon,het_hom[line_list[31].lower()],mute_type,line_list[77],line_list[78],line_list[83],line_list[84],line_list[82],level_value]
                site_Dic_fubiao[line_list[85]].append(site_detail)
                site_li_fubiao.append(site)
    except Exception as e:
            traceback.print_exc(e)
    sort_site = sort_siteDic(site_Dic)
    sort_site_fubiao = sort_siteDic_fubiao(site_Dic_fubiao)
    return disease_Dic,sort_site,sort_site_fubiao

def sort_siteDic(Dic):
    sorted_Dic={}
    for key in Dic:
        sorted_Dic[key]=sorted(Dic[key],key=lambda x:x[11])  #### panqi x[11]
    return sorted_Dic

def sort_siteDic_fubiao(Dic):
    sorted_Dic={}
    for key in Dic:
        sorted_Dic[key]=sorted(Dic[key],key=lambda x:x[10])
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
    if level=='致病':level_value=1
    elif level=='可能致病':level_value=2
    elif '临床意义未明1级' in re.sub(' ','',level): level_value=3
    elif '临床意义未明2级' in re.sub(' ','',level): level_value=4
    elif '临床意义未明3级' in re.sub(' ','',level): level_value=5
    elif '临床意义未明4级' in re.sub(' ','',level): level_value=6
    elif '临床意义未明5级' in re.sub(' ','',level): level_value=7
    elif level=='可能良性': level_value=8
    elif level=='良性':level_value=9
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
    'stopgain2':'缺失移码变异',
    'stoploss':'终止缺失变异',
    'stoploss SNV':'终止缺失变异',
    'synonymous SNV':'同义变异',
    'nonsynonymous SNV':'错义变异',
    'SNV':'单核苷酸变异',
    'splicing SNV':'剪切区单核苷酸变异',
    'splicing INDEL':'剪切区变异',
    }
    if mute_Dic.has_key(mute):
        mute_des=mute_Dic[mute]
    else:
        sys.stderr.write(u'请检查filter文件第9列是否为%s'%mute_Dic.keys())
        sys.exit(1)
    return mute_des

def read_GXYdata(disease_file):
    data=disease_file
    disease_Dic=defaultdict(list)
    if os.path.isfile(data) == False:
        sys.stderr.write('%s is not a file' % data)
        sys.exit(1)
    try:
        with open(data, 'rU') as file_handle:
            for line in file_handle:
                line_1=re.sub('%','\%',line)
                line_2=re.sub('_','\_',line_1)
                line_list=[x.strip() for x in line_2.strip().split('\t')]
                if line_list[0]==u'疾病（共15种）（疾病列表中）':continue
                if len(line_list) < 4:
                    sys.stderr.write(u'cs.data 列数不足')
                    sys.exit(1)
                gene=re.sub(' ','',line_list[1])
                subtype_disease=re.sub(' ','',line_list[2])
                key=gene+':'+subtype_disease
                disease_Dic[key]=[line_list[0],line_list[2],line_list[3]]  # 大类+亚型，亚型名，大类名  
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
    for key in sorted(site_Dic,key=lambda x:site_Dic[x][0][11]):
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
	
	

"""
def get_zongshu(site_Dic,disease_Dic):
    #疾病+亚型0，疾病亚型1，疾病大类2，遗传模式3，疾病描述4，临床表现5，致病性6,指导建议7
    #nue0,acid1,chr_site2,exon3,het/hom4,mute_type5,
    #zhibingxing6,literature7,gene_des8,site_des9,advise10,level_value11,true_false12,subtype_disease,genetic_model
    cover_Dic={}
    gene_li=[]
    disease_li=[]
    genetic_model=[]
    true_false_li=[]
    zhibing_li=[]
    het_hom=[]
    site_num=0
    sortDic=site_Dic
    for key in sorted(sortDic,key=lambda x:sortDic[x][0][11]):   ### panqi ???
        for n in sortDic[key]:
            gene_li.append(key)
            site_num += 1
            disease_li.append(n[13].split('，')[0])
            genetic_model.append(n[14])
            true_false_li.append(n[12])
            zhibing_li.append(n[6])
            het_hom.append(n[4])
    zongshu='此次单基因高血压/血钾异常基因检测，包含44个疾病相关基因，筛查15种疾病亚型。'
    true_box=["trues","trues;DM"]          ### panqi 180817
    mode_Dic={'AD':'常染色体显性遗传',
'AR':'常染色体隐性遗传',
'XLD':'X 染色体连锁显性遗传',
'XLR':'X 染色体连锁隐性遗传',
'XL':'X 染色体连锁遗传',
'AD/AR':'常染色体遗传',
'AR/AD':'常染色体遗传',
'AD,AR':'常染色体遗传',
'AR,AD':'常染色体遗传',
'AD，AR':'常染色体遗传',
'AR，AD':'常染色体遗传',
'AD，Smu':'常染色体显性遗传',
'不明':'遗传模式不明'}
    if site_num == 1:
        if mode_Dic.has_key(genetic_model[0]):
            model_new = mode_Dic[genetic_model[0]]
        else:
            sys.stderr.write(u'遗传模式为%s，不在AD AR XLD XLR XL AD/AR AR/AD  AD,AR  AR,AD  AD，AR  AR，AD,不明内'%mode)
            sys.exit(1)

        if true_false_li[0].lower() == 'false' and u'临床意义未明' in zhibing_li[0]:
            if het_hom[0] == u'杂合' and (genetic_model[0]=='AR' or genetic_model[0]=='XLR'):
                zongshu=r'此次检测有1个变异位点，基因名为%s，疾病亚型为%s，遗传模式为%s，位点为%s位点，位点致病性为%s，与数据库信息匹配结果为%s,隐性杂合，需产品手动添加总述'%(gene_li[0],disease_li[0],genetic_model[0],het_hom[0],zhibing_li[0],true_false_li[0])
            else:
                zongshu += r'''检测结果显示您一共发生了1个位点变异，与%s相关。
该病为%s，'''%(disease_li[0],model_new)
                if gene_li[0]=='ARMC5':zongshu+=r'但当有一个等位基因发生突变时，不会引起瘤状组织的产生，只有当另一个等位基因同时发生体细胞突变时，才会促发肿瘤的发生。'
                zongshu += r'''在与其相关的%s基因上检出一个%s的位点变异。
该变异很罕见，目前尚无相关研究报道，因此致病性不明确，但不排除导致患病的可能性，需要结合其他临床信息综合考虑。
如果您已表现出与该病相关的临床症状，则有可能是该变异导致，建议医生做进一步临床确认并采取相应的治疗方案。
如果您并未表现出与该病相关的临床症状，则该变异可能并不致病，建议您无需紧张，定期体检，并保持积极健康的心态！'''%(gene_li[0],het_hom[0])
        elif true_false_li[0].lower() == 'true\_positive' and zhibing_li[0]== u'致病':
            if het_hom[0] == u'杂合'and (genetic_model[0]=='AR' or genetic_model[0]=='XLR'):
                zongshu=r'此次检测有1个变异位点，基因名为%s，疾病亚型为%s，遗传模式为%s，位点为%s位点，位点致病性为%s，与数据库信息匹配结果为%s,隐性杂合，需产品手动添加总述'%(gene_li[0],disease_li[0],genetic_model[0],het_hom[0],zhibing_li[0],true_false_li[0])
            else:
                zongshu += r'''检测结果显示您一共发生了1个位点变异，与%s相关。该病为%s，在与其相关的%s基因上检出一个%s的位点变异。
根据以往数据，该位点致病的可能性很大，如果您已表现出与该病相关的临床症状，则很可能是该变异导致，建议医生做进一步临床确认并采取相应的治疗方案。'''%(disease_li[0],model_new,gene_li[0],het_hom[0])
        elif (true_false_li[0].lower() in true_box) and (zhibing_li[0]==u'致病' or zhibing_li[0]==u'可能致病'):
            if het_hom[0] == u'杂合' and (genetic_model[0]=='AR' or genetic_model[0]=='XLR'):
                zongshu = r'此次检测有1个变异位点，基因名为%s，疾病亚型为%s，遗传模式为%s，位点为%s位点，位点致病性为%s，与数据库信息匹配结果为%s,隐性杂合，需产品手动添加总述'%(gene_li[0],disease_li[0],genetic_model[0],het_hom[0],zhibing_li[0],true_false_li[0])
            else:
                zongshu += r'''检测结果显示您一共发生了1个位点变异，与%s相关。该病为%s，在与其相关的%s基因上检出一个%s的位点变异。
根据研究报道，该位点致病的可能性较大，如果您已表现出与该病相关的临床症状，则很可能是该变异导致，建议医生做进一步临床确认并采取相应的治疗方案。'''%(disease_li[0],model_new,gene_li[0],het_hom[0])
        elif (true_false_li[0].lower() in true_box) and (r'临床意义未明' in zhibing_li[0]) :
            if het_hom[0] == u'杂合' and (genetic_model[0]=='AR' or genetic_model[0]=='XLR'):
                zongshu = r'此次检测有1个变异位点，基因名为%s，疾病亚型为%s，遗传模式为%s，位点为%s位点，位点致病性为%s，与数据库信息匹配结果为%s,隐性杂合，需产品手动添加总述'%(gene_li[0],disease_li[0],genetic_model[0],het_hom[0],zhibing_li[0],true_false_li[0])
            else:
                zongshu += r'''检测结果显示您一共发生了1个位点变异，与%s相关。该病为%s，'''%(disease_li[0],model_new)
                if gene_li[0]=='ARMC5':zongshu+=r'但当有一个等位基因发生突变时，不会引起瘤状组织的产生，只有当另一个等位基因同时发生体细胞突变时，才会促发肿瘤的发生。'
                zongshu += r'''在与其相关的%s基因上检出一个%s的位点变异。
该变异很罕见，目前研究较少，因此致病性不明确，但不排除导致患病的可能性，需要结合其他临床信息综合考虑。
如果您已表现出与该病相关的临床症状，则有可能是该变异导致，建议医生做进一步临床确认并采取相应的治疗方案。
如果您并未表现出与该病相关的临床症状，则该变异可能并不致病，建议您无需紧张，定期体检，并保持积极健康的心态！'''%(gene_li[0],het_hom[0])
        elif (true_false_li[0].lower() in true_box)  and  (zhibing_li[0]==u'良性' or zhibing_li[0]==u'可能良性'):
            if het_hom[0] == u'杂合' and (genetic_model[0]=='AR' or genetic_model[0]=='XLR'):
                zongshu = r'此次检测有1个变异位点，基因名为%s，疾病亚型为%s，遗传模式为%s，位点为%s位点，位点致病性为%s，与数据库信息匹配结果为%s,隐性杂合，需产品手动添加总述'%(gene_li[0],disease_li[0],genetic_model[0],het_hom[0],zhibing_li[0],true_false_li[0])
            else:
                zongshu += r'''检测结果显示您一共发生了1个位点变异，与%s相关。该病为%s，在与其相关的%s基因上检出一个杂合的位点变异。
但是根据研究报道，该位点为%s变异，致病的可能性不大，建议医生结合临床信息综合判断。'''%(disease_li[0],model_new,gene_li[0],het_hom[0],zhibing_li[0])
        else:
            zongshu =r'此次检测有1个变异位点，基因名为%s，疾病亚型为%s，遗传模式为%s，位点为%s位点，位点致病性为%s，与数据库信息匹配结果为%s,需产品手动添加总述'%(gene_li[0],disease_li[0],genetic_model[0],het_hom[0],zhibing_li[0],true_false_li[0])
    elif site_num == 2:
        if (u'临床意义未明' in zhibing_li[0]) and (u'临床意义未明' in zhibing_li[1]):
            if len(set(gene_li))>1 and len(set(disease_li))>1:
                if het_hom[0] == u'杂合' and het_hom[1]==u'杂合' and (genetic_model[0]=='AR' or genetic_model[0]=='XLR') and (genetic_model[1]=='AR' or genetic_model[1]=='XLR'):
                    zongshu = r'此次检测有2个变异位点，不同基因，不同疾病，隐性杂合，需产品手动添加总述'
                else:
                    zongshu +=r'检测结果显示您一共发生了2个位点变异，与两种疾病相关，'
                    if mode_Dic[genetic_model[0]]==mode_Dic[genetic_model[1]]:
                        zongshu += r'分别为%s和%s。这两种病均为%s。'%(disease_li[0],disease_li[1],mode_Dic[genetic_model[0]])
                    else:
                        zongshu += r'分别为%s（%s）和%s（%s）。'%(disease_li[0],mode_Dic[genetic_model[0]],disease_li[1],mode_Dic[genetic_model[1]])
                    ye_tag=''
                    if len(set(het_hom))==1:
                           ye_tag='也'
                    zongshu += r'其中，在与%s相关的%s基因上检出一个%s的位点变异。在与%s相关的%s基因上%s检出一个%s的位点变异。'%(disease_li[0],gene_li[0],het_hom[0],disease_li[1],gene_li[1],ye_tag,het_hom[1])
                    zongshu += r'这两个变异均很罕见，目前致病性不明确，但不排除导致患病的可能性，需要结合其他临床信息综合考虑。如果您已表现出与上述某种疾病相关的临床症状，则有可能是相关变异导致。建议医生做进一步临床确认并采取相应的治疗方案。如果您并未表现出与上述两种疾病相关的临床症状，则这两个变异可能并不致病。'
            elif len(set(gene_li))>1 and len(set(disease_li))==1:
                if (genetic_model[0]=='AR' or genetic_model[0]=='XLR') and (genetic_model[1]=='AR' or genetic_model[1]=='XLR'):
                    zongshu =r'有2个变异位点，为不同基因，同一疾病，均为隐性杂合或隐性纯合，需产品手动添加总述'
                else:
                    zongshu += r'检测结果显示您一共发生了2个位点变异，与%s相关。该病为%s，在与其相关的%s基因和%s基因上'%(disease_li[0],mode_Dic[genetic_model[0]],gene_li[0],gene_li[1])
                    if len(set(het_hom))==1:
                        zongshu += r'检出两个%s的位点变异。'%het_hom[0]
                    else:
                        zongshu += r'分别检出一个纯合和一个杂合的位点变异。'
                    zongshu += r'这两个变异均很罕见，目前致病性不明确，但不排除导致患病的可能性，需要结合其他临床信息综合考虑。如果您已表现出与该病相关的临床症状，则有可能是上述某一变异导致，建议医生做进一步临床确认并采取相应的治疗方案。如果您并未表现出与该病相关的临床症状，则这两个变异可能并不致病，建议您无需紧张，定期体检，并保持积极健康的心态！'
            elif len(set(gene_li))==1 and len(set(disease_li))==1:
                if het_hom[0] == u'杂合' and het_hom[1]==u'杂合' and (genetic_model[0]=='AR' or genetic_model[0]=='XLR') and (genetic_model[1]=='AR' or genetic_model[1]=='XLR'):
                    zongshu += r'''检测结果显示您一共发生了2个位点变异，与%s相关。
该病为%s，理论上必须在两条等位染色体上同时出现致病性变异才有可能致病（纯合或复合杂合）。此次在与其相关的%s基因上检出两个杂合的位点变异。
这两个变异均很罕见，目前致病性不明确，但不排除导致患病的可能性，需要结合其他临床信息综合考虑。
且这两个变异是位于同一染色体还是分别位于两条等位染色体上需通过家系才能判定，建议联系家人一起进行检测。'''%(disease_li[0],mode_Dic[genetic_model[0]],gene_li[0])
                else:
                    zongshu += r'检测结果显示您一共发生了2个位点变异，与%s相关。该病为%s，在与其相关的%s基因上检出'%(disease_li[0],mode_Dic[genetic_model[0]],gene_li[0])
                    if len(set(het_hom))==1:
                        zongshu += r'两个%s的位点变异。'%het_hom[0]
                    else:
                        zongshu += r'一个%s和一个%s的位点变异。'%(het_hom[0],het_hom[1])
                    zongshu += r'''
这两个变异均很罕见，目前致病性不明确，但不排除导致患病的可能性，需要结合其他临床信息综合考虑。
如果您已表现出与该病相关的临床症状，则有可能是上述某一变异导致，建议医生做进一步临床确认并采取相应的治疗方案。
如果您并未表现出与该病相关的临床症状，则这两个变异可能并不致病，建议您无需紧张，定期体检，并保持积极健康的心态！'''
            else:
                zongshu = r'有2个变异位点，均为临床意义未明，需产品手动添加总述'
        else:
            zongshu = r'有2个变异位点，致病性分别为%s和%s，需产品手动添加总述'%(zhibing_li[0],zhibing_li[1])
    else:
        zongshu = r'位点数为%s,需产品手动添加总述'%str(site_num)
    if '需产品手动添加' in zongshu:
        print '需产品手动添加报告总述'
    else:
        zongshu += r'具体变异信息详见下表，请结合临床分析，如想进一步了解家族遗传史，建议联系家人一起进行检测。'
    cover_Dic['zongshu']=zongshu
    return cover_Dic
"""
	
def get_disease_table(Dic,siteDic):
    #疾病+亚型0，疾病亚型1，疾病大类2，遗传模式3，疾病描述4，临床表现5，致病性6,指导建议7
    #nue0,acid1,chr_site2,exon3,het/hom4,mute_type5,
    #zhibingxing6,literature7,gene_des8,site_des9,advise10,level_value11,true_false12
    cover_Dic={}
    ADAR=[]
    table=''
    sortDic=siteDic
    for key in sorted(sortDic,key=lambda x:sortDic[x][0][11]):
        if len(Dic[key])==1:
            ADAR.append(Dic[key][0][3])
            table += r'''
%s & \makecell[c]{%s} & %s & %s \\ \hline'''%(key,Dic[key][0][0],Dic[key][0][3],Dic[key][0][5])
        else:
            n=0
            for li in Dic[key]:
                ADAR.append(li[3])
                n += 1
                if n==1:
                    table += r'''
\multirow{%s}{*}{%s} & \makecell[c]{%s} & %s & %s \\ \cline{2-4}'''%(str(len(Dic[key])),key,li[0],li[3],li[5])
                elif n==len(Dic[key]):
                    table += r'''
& \makecell[c]{%s} & %s & %s \\ \hline'''%(li[0],li[3],li[5])
                else:
                    table += r''' 
& \makecell[c]{%s} & %s & %s \\ \cline{2-4}'''%(li[0],li[3],li[5])
    genetic_model = get_model(set(ADAR))
    cover_Dic['result_table']=table
    cover_Dic['genetic_model']=genetic_model
    return cover_Dic

def prepare_lite(li):
    li_new=[]
    for i in li:
        if '|' in i:
            arr=i.split('|')
	    for k in arr:
		if k=='' or k=='-' or k=='.':
                	pass
            	else:
			li_new.append(k)
        else:
            li_new.append(i)
    return li_new

def get_model(li):
    genetic_model=[]
    mode_Dic={'AD':'AD：常染色体显性遗传',
    'AR':'AR：常染色体隐性遗传',
    'XLD':'XLD：X 染色体连锁显性遗传',
    'XLR':'XLR：X 染色体连锁隐性遗传',
    'XL':'XL：X 染色体连锁遗传',
    'AD/AR':'AD/AR：常染色体遗传',
    'AR/AD':'AR/AD：常染色体遗传',
    'AD,AR':'AD/AR：常染色体遗传',
    'AR,AD':'AR/AD：常染色体遗传',
    'AD，AR':'AD/AR：常染色体遗传',
    'AR，AD':'AR/AD：常染色体遗传',
    'AD，Smu':'AD：常染色体显性遗传；Smu：体细胞突变',
    '不明':'遗传模式不明'
    }
    ad_f=0
    for i in li:
        mode=re.sub(' ','',i)
	if mode=='AD':
	    if ad_f==0:
		ad_f=1
		genetic_model.append( mode_Dic[mode])
	elif mode == 'AD，Smu':
	    if ad_f==0:
		ad_f=1
		genetic_model.append( mode_Dic[mode])
	    else:
		genetic_model.append('Smu：体细胞突变')
        elif mode_Dic.has_key(mode):
            genetic_model.append( mode_Dic[mode])
        else:
            sys.stderr.write(u'遗传模式为%s，不在AD AR XLD XLR XL AD/AR AR/AD AD,AR AR,AD,不明 内'%mode)
            sys.exit(1)
    genetic='；'.join(genetic_model)
    return genetic

def get_site_table(Dic):
    table=''
    pep=''
    cover_Dic={}
    num=0
    sortDic=Dic
    for key in sorted(sortDic,key=lambda x:sortDic[x][0][11]):
        if len(sortDic[key])==1:
            color=get_color(sortDic[key][0][6])
            pep = sortDic[key][0][1]
            if "fs" in sortDic[key][0][1]:
                pep = "\makecell*[c]{%s}" % sortDic[key][0][1].replace("fs","\\\\fs") #wph add 190612
            table += r'''
\color{MyFontGray}
%s  & %s & %s & %s &  \makecell*[c]{%s} & %s & %s & \textcolor{%s}{%s} \\ \hline
'''%(key,sortDic[key][0][0],pep,sortDic[key][0][2],sortDic[key][0][3],sortDic[key][0][4],sortDic[key][0][5],color,sortDic[key][0][6])
        else:
            num=get_multline(sortDic[key])
            m=0
            for line in sortDic[key]:
                m += 1
                color=get_color(line[6])
                pep=line[1]
                if "fs" in line[1]:
                    pep= "\makecell*[c]{%s}" % line[1].replace("fs","\\\\fs")  #wph add 190612
                if m==1:
                    table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s} & %s & %s & %s &  \makecell*[c]{%s} & %s & %s & \textcolor{%s}{%s} \\ \cline{2-8}'''%(str(num),key,line[0],pep,line[2],line[3],line[4],line[5],color,line[6])
                elif m==len(Dic[key]):
                    table += r'''
\color{MyFontGray}
& %s & %s & %s &  \makecell*[c]{%s} &  %s & %s & \textcolor{%s}{%s} \\ \hline'''%(line[0],pep,line[2],line[3],line[4],line[5],color,line[6])
                else:
                    table += r'''
\color{MyFontGray}
& %s & %s & %s &  \makecell*[c]{%s} &  %s & %s & \textcolor{%s}{%s} \\ \cline{2-8}'''%(line[0],pep,line[2],line[3],line[4],line[5],color,line[6])
    cover_Dic['site_table']=table
    return cover_Dic


def get_fubiao_table(siteDic):
    cover_Dic={}
    table=''
    pep=''
    sortDic=siteDic
    biaoshi=0
    Dic_wenxan=[]
    for key in sorted(sortDic,key=lambda x:sortDic[x][0][10]):
    #for key in  sortDic:
        biaoshi=1
        m1=len(sortDic[key])
        if len(sortDic[key])==1:
            pep = sortDic[key][0][1]
            if "fs" in sortDic[key][0][1]:
                pep = "\makecell*[c]{%s}" % sortDic[key][0][1].replace("fs","\\\\fs")
            table += r'''
\color{MyFontGray}
%s  & %s & %s &  \makecell*[c]{%s} & %s & %s & %s & \makecell*[c]{%s} & %s  & %s  \\ \hline
'''%(key,sortDic[key][0][0],pep,sortDic[key][0][2],sortDic[key][0][3],sortDic[key][0][4],sortDic[key][0][5],sortDic[key][0][6].replace("临床意义","临床意义\\\\"),sortDic[key][0][7],sortDic[key][0][8])
            if sortDic[key][0][9] not in Dic_wenxan:
                Dic_wenxan.append(sortDic[key][0][9]) 
        else:
            m=0;
            num=len(sortDic[key])
            for line in sortDic[key]:
                m += 1
                pep=line[1]
                if "fs" in line[1]:
                    pep= "\makecell*[c]{%s}" % line[1].replace("fs","\\\\fs")
                if  m==1:
                    table += r'''
\color{MyFontGray}
\multirow{%s}{*}{%s}  & %s & %s &  \makecell*[c]{%s} & %s & %s & %s & \makecell*[c]{%s} & %s  & %s  \\  \cline{2-10}
'''%(str(num),key,line[0],pep,line[2],line[3],line[4],line[5],line[6].replace("临床意义","临床意义\\\\"),line[7],line[8])
                elif m==num:
                    table += r'''
\color{MyFontGray}
  & %s & %s &  \makecell*[c]{%s} & %s & %s & %s & \makecell*[c]{%s} & %s  & %s  \\  \hline
'''%(line[0],pep,line[2],line[3],line[4],line[5],line[6].replace("临床意义","临床意义\\\\"),line[7],line[8])
                else:
                    table += r'''
\color{MyFontGray}
  & %s & %s &  \makecell*[c]{%s} & %s & %s & %s & \makecell*[c]{%s} & %s  & %s  \\  \cline{2-10}
'''%(line[0],pep,line[2],line[3],line[4],line[5],line[6].replace("临床意义","临床意义\\\\"),line[7],line[8])
                if line[9] not in Dic_wenxan:
                    Dic_wenxan.append(line[9])    
    biaotou=r'''\newpage
\includegraphics[height=0.4cm]{XKmini.pdf} \zihao{4}{\color{DarkBlue} \sym{ 附表1：}}\par
\par \vspace*{6mm}
\begin{spacing}{1.3}
\color{MyFontGray} \zihao{5}
下表为无ACMG证据、证据矛盾（偏致病类证据与偏良性证据共存）或偏良性证据的位点（可能与相关疾病出现的表型有关，具体情况请结合临床综合分析）
\par \vspace*{6mm}
\arrayrulecolor{MyFontGray}
\renewcommand\arraystretch{1.4}
\fontsize{7}{8.4}\selectfont \color{MyFontGray}
\tablefirsthead{
\rowcolor{DarkBlue}  
\color{white}{\sym{基因}} & \color{white}{\sym{ 核苷酸变异}} & \color{white}{\sym{氨基酸变异}}  &  \makecell[c]{\color{white} \sym{转录本}\\ \sym{外显子编号}}  & \color{white}{\sym{变异状态}} & \color{white}{\sym{变异类型}} &  \makecell[c]{\color{white} \sym{ACMG} \\\sym{条目}}  &  \color{white}{\sym{位点致病性}} & \color{white}{\sym{相关疾病}} & \color{white}{\sym{遗传模式}} \\}
\tablehead{}
\tabletail{\hline}
\tablelasttail{}
\begin{supertabular}{| m{9mm}<{\centering} | m{15mm}<{\centering} | m{13mm}<{\centering} | m{15mm}<{\centering} | m{6mm}<{\centering} | m{10mm}<{\centering} | m{10mm}<{\centering} |m{13mm}<{\centering}| m{31mm}<{\centering} | m{6mm}<{\centering} | }
'''
    biaowei=r'''\end{supertabular}
%\par \vspace*{3mm}
%\zihao{6}
%注：关于频率，“.”表示在对应数据库中未收录该位点
\end{spacing}'''
    fubiao_tablea=''
    fubiao_wenxiana=''
    if biaoshi==1:
        fubiao_tablea=biaotou+table+biaowei
        fubiao_wenxiana=get_wenxianfu(Dic_wenxan)
    cover_Dic['fubiao_table']=fubiao_tablea
    cover_Dic['fubiao_wenxian']=fubiao_wenxiana
    return cover_Dic


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
    #nue0,acid1,chr_site2,exon3,het/hom4,mute_type5,
    #zhibingxing6,literature7,gene_des8,site_des9,advise10,level_value11
    site_sortDic=site
    disease_Dic=disease
    cover_Dic={}
    jiexi=''
    table_li=[]
    jiexi_beizhu=r'''
\par \\   %vspace*{2mm}
\color{DarkBlue}{\zihao{6} 备注}\\
\zihao{-6} \color{MyFontGray} 预测软件包括 SIFT、PolyPhen2、M-CAP和REVEL, 用于预测错义变异是否会导致蛋白结构和功能发生改变, 准确率大致在 65\%-80\%。需要说明的是, 各软件提供的只是基于算法的预测结果, 可以为位点解读提供参考, 但不足以作为致病性判定的唯一标准, 还需要结合其他信息综合判断。\\ \\
以上位点变异为基因检测的客观结果，但变异是否真正致病需要由医生结合临床信息综合判断。如果您想进一步了解基因变异的信息或家族遗传史 ，建议联系家人一起进行检测。'''
    jiexi_des_head=r'''
\vspace*{-4mm}
\begin{spacing}{1.5}
\color{MyFontGray} \zihao{-5}
\rowcolors{1}{}{LightBlue}
\arrayrulecolor{DarkBlue}
\renewcommand\arraystretch{1.5}
\begin{longtable}{C{3.5cm} L{12.7cm}}'''
    jiexi_des_foot=r'''
\hline
\end{longtable}
\end{spacing}'''
    for key in sorted(site_sortDic,key=lambda x:site_sortDic[x][0][11]):
        jiexi=''
        disease_li,disease_des,site_li,site_des,level=prepare_jiexi(site_sortDic[key],disease_Dic[key])
        #site_li,site_des,level=prepare_jiexi(site_sortDic[key])
        if len(site_li)==1:
            #if disease[key][0][2]==disease[key][0][1]:type_sub=disease[key][0][2]
            #else:type_sub=disease[key][0][2]+'（'+disease[key][0][1]+'）'
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
\end{spacing} '''%(key,site_li[0],'、'.join(disease_li).replace("\\",""),level[0])
            jiexi += jiexi_des_head
            jiexi += r'''
变异解析 & %s \\'''%site_des[0]
            jiexi += r'''
基因描述 & %s \\'''%site[key][0][8]
#            jiexi += r'''
#\makecell[{{p{3.4cm}}}]{%s} & %s\\'''%(type_sub,disease[key][0][4])
            jiexi += get_disease_row(disease_li,disease_des)
            jiexi += jiexi_des_foot
            jiexi += get_wenxian(site[key])
            jiexi += jiexi_beizhu
            jiexi += get_advise(disease[key])
            table_li.append(jiexi)
        else:
            jiexi += get_head(key,site_li,disease[key],level)
            jiexi += jiexi_des_head
            for index,item in enumerate(site_des):
                jiexi +=r'''
变异解析%s & %s \\'''%(str(index+1),item)
            jiexi += r'''
基因描述 & %s \\'''%site[key][0][8]
            """for i in disease[key]:
                if i[2]==i[1]:type_sub=i[2]
                else:type_sub=i[2]+'（'+i[1]+'）'
                jiexi += r'''
\makecell[{{p{3.4cm}}}]{%s} & %s\\'''%(i[0],i[4])"""
            jiexi += get_disease_row(disease_li,disease_des)
            jiexi += jiexi_des_foot
            jiexi += get_wenxian(site[key])
            jiexi += jiexi_beizhu
            jiexi += get_advise(disease[key])
            table_li.append(jiexi)
    table=r'\newpage '.join(table_li)
    cover_Dic['jiexi_table']=table
    return cover_Dic

def get_advise(li):
    advise=r'\par \vspace*{5mm}'
    for i in li:
        if i[6]==u'致病' or i[6]==u'可能致病':
            advise += r'''
\includegraphics[height=0.4cm]{XKmini.pdf} \zihao{4}{\color{DarkBlue}\sye \sym{ %s治疗指导建议}}\par
\fontsize{10.5}{13.65}\selectfont{\color{MyFontGray} %s}\par'''%(i[2],i[7])
    return advise

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
%\par \vspace*{5mm}'''
    else:
        lite =r'''
\par \vspace*{3mm}
\zihao{6}{\color{DarkBlue} 参考文献}\\ \vspace*{-6mm}
\begin{spacing}{1.5}
\zihao{-6} \color{MyFontGray}
\setdefaultleftmargin{2em}{}{}{}{.5em}{.5em}
\begin{compactenum}'''
        for i in wenxian_new:
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
            if i=='' or i=='-' or i=='.':
                pass
            else:
                wenxian.append(i)
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


def get_head(key,site,disease,level):
    head=''
    disease_new=[]
    for i in disease:
        if i[2]==i[1]:type_sub=i[2]
        else:type_sub=i[2]+'（'+i[1]+'）'
        disease_new.append(type_sub)
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
        for index,item in enumerate(site):
            head +=r'''
风险位点%s: & %s—— %s, %s\\'''%(str(index+1),key,item,level[index])
        head +=r'''
相关疾病: & %s\\
\end{tabular}}
\end{spacing} '''%('、'.join(disease_new))
    return head

def get_disease_row(disease,disease_des):
    table=''
    if len(disease) != len(disease_des):
        sys.stderr.write(u'请检查filter是否每个位点都包含疾病名称、疾病描述')
        sys.exit(1)
    else:
        for index,item in enumerate(disease):
            table += r'''
%s & %s \\'''%(item.replace("\\",""),disease_des[index])
    return table

def prepare_jiexi(site_li,disease_li): #wph changed 190612
    #nue,acid,chr_site,exon,het/hom/,mute_type,zhibingxing,literature,gene_des,site_des
    disease=[]
    disease_des=[]
    site=[]
    site_des=[]
    level=[]
    for line in disease_li:
        if line[0] not in disease:
            disease.append(line[0])
            disease_des.append(line[4])
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
    if module.strip().upper() == 'GXY':
        cover_Dic['title'] = 'Part'+part+' — 检测结果'
    else:
        cover_Dic['title'] = 'Part'+part+' — 单基因高血压/血钾异常基因检测结果'
    out_GXY=open(outDir+'/'+sampleid+'_Cardio_GXY.tex','w')
    if os.path.isfile(yinyang_file) == False:
        sys.stderr.write('%s is not a file' % yinyang_file)
        sys.exit(1)
    try:
        yinyang_handle = open(yinyang_file,'r')
        yinyang = yinyang_handle.readlines()[0].strip()
        yinyang_handle.close()
        if yinyang == '0':
            id_filter,id_filterfu = get_filter(filterFile,sampleid)
            with open('/PUBLIC/pipline/script/Report/XKFW/V3_180824/tempt/Cardio_GXYyin_temple.tex','r') as file_handle:
                GXY_temp = Template(file_handle.read())
                diseaseDic,siteDic,siteDic_fu = prepare_allDic(id_filter,disease_file,id_filterfu)
                fubiao_table = get_fubiao_table(siteDic_fu)
                cover_Dic.update(fubiao_table)
                num=len(id_filterfu)
                if num==0:
                    cover_Dic['yin']=r'''此次单基因高血压/血钾异常基因检测，包含44个疾病相关基因，筛查15种疾病亚型。检测结果显示您并未发生明确致病变异的位点，因此，您的疾病症状并不是由这些已知基因的致病变异引起的，可能是目前研究尚未涉及的基因变异，也可能是由于其他因素，如不良生活习惯、环境、其他疾病等引起。建议您注意生活细节，积极改变不良生活方式，并配合医生治疗。 

以上检测结果仅基于当前的科学研究，随研究文献的更新，疾病-基因的对应关系可能会发生变化，因此我们的报告具有一定的时效性 。
'''             
                else:
                    cover_Dic['yin']=r'''此次单基因高血压/血钾异常基因检测，包含44个疾病相关基因，筛查15种疾病亚型。检测结果显示您并未发生高风险的位点变异；一些变异位点由于相关研究较少、ACMG证据条目不足，基于现有研究判断其致病性等级较低，但不排除导致疾病的可能性，这类位点的致病性等级及相关证据条目展示在附表中，供您和临床医生参考。

若检出位点的相关疾病与您现有的临床表型不符，那么您的相关临床症状可能与其他未知基因变异、疾病的不完全外显性、不良生活习惯或环境等因素相关。建议您注意生活细节，积极改变不良生活方式，或去医院进行相关疾病的进一步检测，听从医生的临床建议，定期体检，并采取积极的生活干预。

若您未发现任何相关疾病表型，则这类 ACMG 致病类证据较少的位点可能并不致病。 

以上检测结果仅基于当前的科学研究，随研究文献的更新，位点致病性情况可能会发生变化，属于正常现象。如果您想进一步了解基因变异的信息或家族遗传史，建议联系家人一起进行检测。
'''
                cover_Dic['num']=num
                #module_li=sorted([x.strip().upper() for x in module.strip().split('+')])
                #gxy_module=['GXY','DRUGB']
                #module_gxy=sorted(['GXY','DRUGB'])
                #module_gxy2=sorted(['GXY','DRUGB','CS'])
                #if module_li == module_gxy: cover_Dic['yin']=r'''
#\fontsize{10.5}{13.65}\selectfont{\color{MyFontGray} 为了使您更有效的服用降压药，我们还为您检测了 40 种心血管常用药物的敏感基因位点，涵盖降压药、降脂药、降糖药和抗凝药，并给出了具体用药指导,请您详见 Part3 部分报告。}\\ \\'''
                #elif module_li == module_gxy2: cover_Dic['yin']=r'''
#\fontsize{10.5}{13.65}\selectfont{\color{MyFontGray} 为了使您更有效的服用降压药，我们还为您检测了 40 种心血管常用药物的敏感基因位点，涵盖降压药、降脂药、降糖药和抗凝药，并给出了具体用药指导,请您详见 Part4 部分报告。}\\ \\'''
                #else: cover_Dic['yin']=''
                report = GXY_temp.safe_substitute(cover_Dic)
                print >> out_GXY,report
        elif yinyang == '1':
            id_filter,id_filterfu = get_filter(filterFile,sampleid)
            with open('/PUBLIC/pipline/script/Report/XKFW/V3_180824/tempt/Cardio_GXYyang_temple.tex','r') as file_handle:
                GXY_temp = Template(file_handle.read())
                diseaseDic,siteDic,siteDic_fu = prepare_allDic(id_filter,disease_file,id_filterfu)
                #zongshu = get_zongshu(siteDic,diseaseDic)
                zongshu = get_zongshu(siteDic)
                disease_table = get_disease_table(diseaseDic,siteDic)
                site_table = get_site_table(siteDic)
                fubiao_table = get_fubiao_table(siteDic_fu)
                jiexi_table = get_jiexi_table(diseaseDic,siteDic)
                cover_Dic.update(zongshu)
                cover_Dic.update(disease_table)
                cover_Dic.update(site_table) 
                cover_Dic.update(fubiao_table)
                cover_Dic.update(jiexi_table)
                report = GXY_temp.safe_substitute(cover_Dic)
                print >> out_GXY,report
    except Exception as e:
        traceback.print_exc(e)
    finally:
        out_GXY.close()

if __name__ == '__main__':
    main()
