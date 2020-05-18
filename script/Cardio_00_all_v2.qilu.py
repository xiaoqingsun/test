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
import traceback
from datetime import datetime
from xlrd import xldate_as_tuple

def argInit():
    parser = argparse.ArgumentParser(
        description="report_create program.")
    parser.add_argument(
        '-cfg', help='file cfg.txt', required=True)
    parser.add_argument(
        '-qc', help='file Stat_QC.xls', required=True)
    parser.add_argument(
        '-sample', help='sampleid', required=True)
    parser.add_argument(
        '-feature', help='feature', required=True)
    parser.add_argument(
        '-cs', help='猝死的阴阳性判定文件cs.yinYang.txt', required=True)
    parser.add_argument(
        '-gxy', help='高血压的阴阳性判定文件gxy.yinYang.txt', required=True)
    parser.add_argument(
        '-module', help='CS+GXY+Drug or CS or GXY or ^', required=True)
    parser.add_argument(
        '-out', help='project name', default=os.getcwd())                   ##########缺失路径
    argv = vars(parser.parse_args())
    global QCfile, cfgFile, sampleid, Module,outDir,module,cs_file,gxy_file,feature_file
    QCfile = argv['qc'].strip()
    cfgFile = argv['cfg'].strip()
    sampleid=argv['sample'].strip()
    Module = argv['module'].strip()
    cs_file = argv['cs'].strip()
    gxy_file = argv['gxy'].strip()
    feature_file = argv['feature'].strip()
    outDir = os.path.abspath(argv['out'].strip())
    module = module_zuhe(Module)

def module_zuhe(Module):
    module_1=re.sub('Druga','Drug',Module)
    module_2=re.sub('Drugb','Drug',module_1)
    module='+'.join(set(module_2.split('+')))
    return module

def get_Title_Letter(module,Module,title_file):
    #title_file='/PUBLIC/work/sunxiaoqing/python_word/QiLuHosputal/cover.txt'
    print module
    print Module
    cover_Dic = defaultdict(list)
    module_li=sorted([x.strip().upper() for x in module.strip().split('+')])
    Module_li=sorted([x.strip().upper() for x in Module.strip().split('+')])
    if os.path.isfile(title_file) == False:
        sys.stderr.write('%s is not a file' % title_file)
        exit
    try:
        with open(title_file, 'rU') as file_handle:
            for line in file_handle:
                line_list = line.strip().split('\t')
                arr = sorted([x.strip().upper() for x in line_list[0].split('+')])
                if module_li == arr or Module_li == arr:
                    cover_Dic['Title']=line_list[1]
                    cover_Dic['Letter']=re.sub(r'%',r'\%',line_list[2])
                    cover_Dic['Part']=line_list[3]
                    cover_Dic['header']=line_list[4]
                    cover_Dic['Fuwuliucheng']=line_list[5]
                    break
                else:
                    pass
                if cover_Dic['Title']=='' or cover_Dic['Letter']=='':
                    sys.stderr.write('please check model or file %s' % title_file)
                    exit
    except Exception, e:
        traceback.print_exc(e)
    return cover_Dic

def get_sampleInf(sampleid,cfgFile,feature_file):
    sample_inf_list=[]
    hospital=''
    cover_Dic = defaultdict(list)
    times_now = a=time.strftime("%Y-%m-%d", time.localtime())
    head=['id','name','sex','age','date1','chenghu','feature','date2']
    if os.path.isfile(cfgFile) == False:
        print('%s is not a file' % cfgFile)
        return cover_Dic
    try:
        with open(cfgFile, 'rU') as file_handle:
            for line in file_handle:
                line_list=[x.strip() for x in line.strip().split('\t')]
                if line_list[0]==sampleid:
                    line_list[5]=re.sub(r'\(.*?\)','',line_list[5])
                    (age,date,sex)=get_HTB(sampleid) #wph changed 190423 修改性别不匹配的问题
                    hospital=line_list[9]
                    if sex not in line_list[5]:
                        line_list[5] = sex
                    chenghu=judge_chenghu(line_list[5])
                    #(age,date)=get_HTB(sampleid)
                    feature=get_feature(feature_file,sampleid)
                    print (len(line_list[4]))
                    if len(line_list[4])<9:
                        line_list[4]=line_list[4]+'{\\quad}'
                    sample_inf_list=[line_list[0],line_list[4],line_list[5],age,date,chenghu,feature,times_now]  #需要送样日期，暂时用出报告日期代替
                    break
        cover_Dic = {k:v for k,v in zip(head,sample_inf_list)}
    except Exception, e:
        traceback.print_exc(e)
    return cover_Dic,hospital

def judge_chenghu(sex):
    chenghu=''
    try:
        if sex ==r'女':
            chenghu=r'女士'
        elif sex ==r'男':
            chenghu=r'先生'
        else:
            print('性别是%s,check please' % sex)
    except Exception, e:
        traceback.print_exc(e)
    return chenghu

def get_feature(feature_file,sampleid):   ### panqi add 181125
    feature = '无临床信息'
    if os.path.isfile(feature_file) == False:
        print('%s is not a file' % feature_file)
        return 'Error File'
    try:
        file_handle = open(feature_file, 'rU')
        feature_list = file_handle.readlines()
        for line in feature_list:
            feature_arr = [x.strip() for x in line.strip().split('\t')]
            if feature_arr[0] == sampleid and len(feature_arr)>1:
                feature=r''' %s ''' %(feature_arr[1])
            else:
                pass
    except Exception, e:
        traceback.print_exc(e)
    return feature


def get_HTB(sampleid):
    #files=glob.glob(r'/PUBLIC/pipline/database/sheet_hutong/*心康互通表*.xls*')
    files=glob.glob(r'/PUBLIC/pipline/database/sheet_hutong/2019心康送样信息单*.xls*')
    max_num=0
    age=date=''
    sex=''
    age_num = '' #wph changed 181106
    date_num = ''
    for i in files:
        hutong=i.split('/')[-1]
        #st=re.match(r'心康互通表（([0-9]*)）.xlsx',hutong)
        st=re.match(r'2019心康送样信息单-([0-9]*).xlsx',hutong)
        if int(st.group(1))> max_num:
            max_num=int(st.group(1))
    #HTB=r'/PUBLIC/pipline/database/sheet_hutong/心康互通表（%s）.xlsx' %max_num
    HTB=r'/PUBLIC/pipline/database/sheet_hutong/2019心康送样信息单-%s.xlsx' %max_num
    try:
        data = xlrd.open_workbook(HTB)
        table = data.sheet_by_name('下单表')
        nrows = table.nrows
        for rownum in range(1,nrows):
            if table.row(rownum)[2].value == sampleid:
                age=table.row(rownum)[8].value
                sex=table.row(rownum)[7].value
                ctype = table.cell(rownum, 1).ctype
                cell =table.cell_value(rownum, 1)
                date
                if ctype ==3 :
                    date = datetime(*xldate_as_tuple(cell, 0))
                    date = date.strftime('%Y/%m/%d')
                else:
                    date=table.row(rownum)[1].value
        date_num=date.replace("/","-")
                
        if isinstance(age, float) == True: #wph changed 181106
            age_num=unicode(int(age)) + "岁"
        elif not age:
            age_num="不详"
        else:
            age_num=age
        if age_num=='' or date_num=='':
            sys.stderr.write(u"ID不在互通表 或 互通表中没有年龄/送样日期信息")
            exit
    except Exception, e:
        traceback.print_exc(e)
    return age_num,date_num,sex

def get_Shengming(module):
    shengming=''
    cover_Dic= defaultdict(list)
    module_li=[x.strip().upper() for x in module.strip().split('+')]
    if ('DRUG' in module_li) or ('ALL' in module_li):
        shengming=r'''
\zihao{-4}{\color{DarkBlue} \sym \sye {4. 关于药物}} \par
\zihao{-5}{\color{MyFontGray} 本报告的用药建议严格遵守最新的药物基因组学知识库信息。用药建议仅供参考，具体措施请以临床情况和医生
意见为准。由于儿童的肝脏未发育完全，对药物的代谢能力与成人完全不同，故不能直接沿用成人标准，请在儿童用药时结合医生意见谨慎参考。} \par
\zihao{-4}{\color{DarkBlue} \sym \sye {5. 关于隐私保护}} \par
\zihao{-5}{\color{MyFontGray} 在您的基因检测报告中包含了您本人的基因检测数据资料，这些资料属于个人隐私，请您务必妥善保管，以免资
料的泄露有可能对您个人及家庭造成不利影响。我们郑重承诺妥善保管您的相关数据，未经本人授权不得用于其他用途。若因您个人原因发生信息外泄，其后果将由您本人承担。} \par'''
        cover_Dic['Shengming']=shengming
    else:
        shengming=r'''
\zihao{-4}{\color{DarkBlue} \sym \sye {4. 关于隐私保护}} \par
\zihao{-5}{\color{MyFontGray} 在您的基因检测报告中包含了您本人的基因检测数据资料，这些资料属于个人隐私，请您务必妥善保管，以免资
料的泄露有可能对您个人及家庭造成不利影响。我们郑重承诺妥善保管您的相关数据，未经本人授权不得用于其
他用途。若因您个人原因发生信息外泄，其后果将由您本人承担。} \par'''
        cover_Dic['Shengming']=shengming
    return cover_Dic
    pass

def get_QC_Inf(QCfile,sampleid):
    QC_inf_list=[]
    QC_li=[]
    cover_Dic= defaultdict(list)
    head=['cover_one','cover_four','cover_ten','cover_twenty','average_depth']
    if os.path.isfile(QCfile) == False:
        print('%s is not a file' % QCfile)
        return cover_Dic
    try:
        file_handle = open(QCfile, 'rU')
        QC_list = file_handle.readlines()
        for line in QC_list:
            QC_arr = [x.strip() for x in line.strip().split('\t')]
            if QC_arr[0] == sampleid and len(QC_arr) > 40:
                QC_li.append(QC_arr)
            else:
                pass
        if len(QC_li)==1:
            QC_inf = [QC_li[0][2],QC_li[0][34],QC_li[0][33],QC_li[0][32],QC_li[0][3]]
        else:
            sys.stderr.write(u'请检查质控结果是否无或未合并')
            sys.exit(1)
        cover_Dic = {k:v for k,v in zip(head,QC_inf)}
    except Exception, e:
        traceback.print_exc(e)
    finally:
        if 'fileHandle' in dir() or file_handle.closed == False:
            file_handle.close()
    return cover_Dic

def get_Jiance(module):
    module_li=[x.strip().upper() for x in module.strip().split('+')]
    cover_Dic= defaultdict(list)
    Jiance=''
    Jiance_li=[]
    try:
        if len(module_li) == 1:
            if module_li[0] == 'ALL':
                Jiance_li.append(get_GXY('检测项目1'))
                Jiance_li.append(get_CS('检测项目2'))
                Jiance_li.append(get_drug('检测项目3'))
            elif module_li[0] == 'GXY':
                Jiance_li.append(get_GXY('检测项目'))
            elif module_li[0] == 'CS':
                Jiance_li.append(get_CS('检测项目'))
            elif module_li[0] == 'DRUG':
                Jiance_li.append(get_drug('检测项目'))
        else:
            m=0
            if 'GXY' in module_li:
                m += 1
                num = '检测项目' + str(m)
                Jiance_li.append(get_GXY(num))
            if 'CS' in module_li:
                m += 1
                num = '检测项目' + str(m)
                Jiance_li.append(get_CS(num))
            if 'DRUG' in module_li:
                m += 1
                num = '检测项目' + str(m)
                Jiance_li.append(get_drug(num))
        Jiance = '\\newpage '.join(Jiance_li)
        cover_Dic['Jiance'] = Jiance
    except Exception, e:
        traceback.print_exc(e)
    return cover_Dic

def get_CS(part):
    CS_li=[]
    CS=r'''
\begin{center} \zihao{3}{\color{DarkBlue} \sye \sym %s} \end{center}
\par %%\vspace*{6mm}
\fontsize{10.5}{12.6}\selectfont{\color{MyFontGray} 本次遗传性心血管 病基因检测包括 179 个基因的全部外显子区域，共涉及5大类、 36 种遗传性心血管病，
具体疾病名称请参见下表。} \par'''% part
    CS_li.append(read_CS('xinlv'))
    CS_li.append(read_CS('xinji'))
    CS_li.append(read_CS('zhudongmai'))
    CS_li.append(read_CS('xinzang'))
    CS_li.append(read_CS('qita'))
    CS += '\\vspace*{-10mm}'.join(CS_li)
    CS += r'''
\vspace*{-3mm}
\fontsize{9}{11.7}\selectfont{\color{MyFontGray}注：以上检测疾病中无“*”标识的疾病与基因间的相关性均来自于国内外指南或专家共识；带“*”标识的疾病与基因间的相关
性仅来自于当前最新的科学研究，尚未在指南和共识中明确指出。}'''
    return CS

def read_CS(disease_type):
    CS_Dic=defaultdict(list)
    head_Dic={'xinlv':'遗传性心律失常（59个）',
    'xinji':'遗传性心肌病（101个）',
    'zhudongmai':'遗传性主动脉病（12个）',
    'xinzang':'先天性心脏病（24个）',
    'qita':'其他心脏相关疾病（12个）'
    }
    CS_file='/PUBLIC/pipline/script/Report/XKFW/V3_180824/database/cs_disease_'+ disease_type+'.txt'
    CSxiang=r'''
\begin{spacing}{1.9}
\renewcommand\arraystretch{1.4}
\color{MyFontGray} \zihao{-5}
\arrayrulecolor{DarkBlue}
\rowcolors{2}{}{LightBlue}
\begin{longtable}{L{9cm} L{6.9cm} }
\rowcolor{DarkBlue}
\multicolumn{2}{l}{\zihao{-4} \color{white} \sym \sye{%s}}\\
\zihao{-4}{\color{MyFontGray} \sym{检测疾病}} & \zihao{-4}{\color{MyFontGray} \sym{检测基因}} \\ ''' %head_Dic[disease_type]
    if os.path.isfile(CS_file) == False:
        print('%s is not a file' % CS_file)
        return CSxiang
    try:
        with open(CS_file, 'rU') as file_handle:
            for line in file_handle:
                line_list=[x.strip() for x in line.strip().split('\t')]
                CSxiang += r'''
\raisebox{-0.2ex}{\shortstack[l]{%s}} &  %s \\ ''' %(line_list[0],line_list[1])
        CSxiang += r'''
\hline
\end{longtable}
\end{spacing}'''
    except Exception, e:
        traceback.print_exc(e)
    return CSxiang

def get_GXY(part):
    GXY_Dic=defaultdict(list)
    GXY_file='/PUBLIC/pipline/script/Report/XKFW/V3_180824/database/gxy_disease_name3.txt'
    GXYxiang = r'''
\begin{center} \zihao{3}{\color{DarkBlue} \sye \sym %s} \end{center}
\par %%\vspace*{6mm}
\fontsize{10.5}{12.6}\selectfont{\color{MyFontGray} }本检测共筛查15种单基因高血压/血钾异常疾病，包含44个相关基因，具体疾病名称请参见下表。 \par'''%part
    GXYxiang += r'''
\begin{spacing}{1.2}
\renewcommand\arraystretch{2.1}
\color{MyFontGray} \zihao{-5}
\arrayrulecolor{DarkBlue}
\rowcolors{2}{}{LightBlue}
\begin{longtable}{L{8.3cm} L{7.2cm} }
\rowcolor{DarkBlue}
\zihao{-4}{\color{white} \sym{检测疾病}} & \zihao{-4}{\color{white} \sym{致病基因}} \\ '''
    if os.path.isfile(GXY_file) == False:
        print('%s is not a file' % GXY_file)
        return GXYxiang
    try:
        with open(GXY_file, 'rU') as file_handle:
            for line in file_handle:
                line_list=[x.strip() for x in line.strip().split('\t')]
                GXYxiang += r'''
\makecell[{{p{8.2cm}}}]{%s} & \makecell[{{p{7cm}}}]{%s} \\ ''' %(line_list[0],line_list[1])
        GXYxiang += r'''
\hline
\end{longtable}
\end{spacing}'''
    except Exception, e:
        traceback.print_exc(e)
    return GXYxiang

def get_drug(part):
    drug_xiang =r'''
\begin{center} \zihao{3}{\color{DarkBlue} \sye \sym %s} \end{center}
\par %%\vspace*{6mm}
\fontsize{10.5}{12.6}\selectfont{\color{MyFontGray} 本次心血管用药指导基因检测包括44个基因，共涉及13大类、40种心血管常用药，具体用药名称请参照下表。 \par'''%part
    drug_Dic=defaultdict(list)
    count=0
    drug_file='/PUBLIC/pipline/script/Report/XKFW/V3_180824/database/Drug.txt'
    drug_xiang +=r'''
\begin{spacing}{1.3}
\color{MyFontGray}
\arrayrulecolor{DarkBlue}
\renewcommand\arraystretch{1.1}
%\renewcommand{\multirowsetup}{\centering}
\zihao{-5}
\begin{longtable}{L{1.8cm}  L{2.9cm} L{4.5cm} L{6.5cm} }
\rowcolor{DarkBlue}
&   \zihao{-4}{\color{white} \sym{药物种类}}   &   \zihao{-4}{\color{white} \sym{药物名称}}   &  \zihao{-4}{\color{white} \sym{敏感基因}} \\  '''
    if os.path.isfile(drug_file) == False:
        print('%s is not a file' % drug_file)
        return drug_xiang
    try:
        with open(drug_file, 'rU') as file_handle:
            allline=file_handle.readlines()
        num=0
        for line in allline:
            if re.match('^大类',line):continue
            arr=[x.strip() for x in line.strip().split('\t')]
            if len(arr)==4:
                yalei = get_rowcolor(num,arr[1])
                num += 1
                drug_xiang += r'''
\multirow{-%s}{*}{\zihao{5}\sym{%s}} & %s & %s & %s \tabularnewline \hline'''%(arr[0].split(':')[1],arr[0].split(':')[0],yalei,arr[2],arr[3])
            elif len(arr)==3:
                yalei = get_rowcolor(num,arr[0])
                num += 1
                drug_xiang += r'''
 & %s & %s & %s \tabularnewline '''%(yalei,arr[1],arr[2])
            elif len(arr)==2:
                yalei = get_rowcolor(num,'')
                drug_xiang += r'''
 & %s & %s & %s \tabularnewline '''%(yalei,arr[0],arr[1])
        drug_xiang += r'''
\end{longtable}
\end{spacing}'''
    except Exception, e:
        traceback.print_exc(e)
    return drug_xiang

def get_rowcolor(num,yname):
    color=''
    if num%2==1:color=r'\rowcolor{LightBlue}'
    yalei=''
    if yname=='':
        yalei=color
    elif  yname.split(':')[1]=='1':
        yalei=r'%s %s'%(color,yname.split(':')[0])
    elif yname.split(':')[1]!='1':
        yalei=r'%s \multirow{-%s}{2.9cm}{%s}'%(color,yname.split(':')[1],yname.split(':')[0])
    return yalei

def get_insert_module(module,sampleid,gxy,cs,hospital):
    ordertag=0
    CS_Dic={}
    insert=''
    insert2=''
    if '齐鲁' in hospital:
        ordertag=1
    include_gxy=r'\include{%s_Cardio_GXY}'%sampleid
    include_cs=r'\include{%s_Cardio_CS}'%sampleid
    include_drug=r'\include{%s_Cardio_Yaomim}'%sampleid
    include_sanger=r'\include{%s_Cardio_Sanger}'%sampleid
    mingcishuoming=r'''
%%%%%%%%%%%%%%%%%
\newpage
\addcontentsline{toc}{subsection}{ⅲ 名词注释}
\noindent
\begin{center}
\zihao{2}{\color{DarkBlue} \sye \sym 名词注释}
\end{center}
\par \vspace*{10mm}
\begin{spacing}{1.5}
\color{MyFontGray} \zihao{-5}
1. 杂合：指同一位点上的两个等位基因有不相同的基因型。\\
2. 纯合：指同一位点上的两个等位基因有相同的基因型。\\
3. 单碱基变异：\\
（1）同义变异：由于遗传密码子存在简并性，碱基置换后密码子虽然发生改变，但所编码的氨基酸没有改变。同义变异常发生在三联密码子的第3个碱基。\\
（2）错义变异：碱基置换后编码某个氨基酸的密码子变成另一种氨基酸的密码子，从而改变多肽链的氨基酸序列，影响蛋白质的功能。\\
（3）无义变异：碱基置换后使原本编码氨基酸的密码子变成不编码任何氨基酸的终止密码子（UAG、UAA或UGA），使得多肽链的合成提前终止，肽链长度变短而成为无活性的截短蛋白。\\
4. 终止密码子变异：与无义变异相反，碱基替换后使某一终止密码子变成具有氨基酸编码功能的遗传密码子，使本应终止延伸的多肽链合成异常地持续进行。\\
5. 移码变异\\
（1）移码插入变异：由于编码序列中插入的碱基数目为非3的倍数，使得插入点下游的三联密码子组合发生改变，造成变异点下游的全部氨基酸序列发生改变。\\
（2）移码缺失变异：由于编码序列中缺失的碱基数目为非3的倍数，使得缺失点下游的三联密码子组合发生改变，造成变异点下游的全部氨基酸序列发生改变。\\
（3）移码替换变异： 编码序列中的碱基被另一段碱基替换，且替换之后会引起下游的三联密码组合发生改变，造成变异点后的氨基酸序列发生改变。\\
6. 非移码变异\\
（1）非移码缺失变异：由于编码序列中缺失的碱基数目为3的倍数，导致氨基酸链缺失1个或多个氨基酸，但对下游的三联密码子组合无影响。\\
（2）非移码插入变异：由于编码序列中插入的碱基数目为3的倍数，导致氨基酸链增加1个或多个氨基酸，但对下游的三联密码子组合无影响。\\
（3）非移码替换变异：编码序列中的碱基被另一段碱基替换，且替换之后会引起一个或多个氨基酸发生改变，但不会造成下游的三联密码组合发生改变。\\
\end{spacing}'''
    module_li=sorted([x.strip().upper() for x in module.strip().split('+')])
    print module_li[0]
    if len(module_li)==1 and module_li[0]=='ALL':
        if gxy=='1' and cs =='1':
            insert = include_gxy + include_cs + mingcishuoming + +include_drug
            insert2 = include_sanger
            if ordertag ==1: insert = include_gxy + include_cs + mingcishuoming + include_sanger + include_drug
        elif gxy=='1' and cs =='0':
            insert = include_gxy + mingcishuoming + include_cs + include_drug
            insert2 = include_sanger
            if ordertag ==1: insert = include_gxy + mingcishuoming + include_cs + include_sanger + include_drug
        elif gxy=='0' and cs =='1':
            insert = include_gxy + include_cs + mingcishuoming + include_drug
            insert2 = include_sanger
            if ordertag ==1: insert = include_gxy + include_cs + mingcishuoming + include_sanger +include_drug
        elif gxy=='0' and cs =='0':
            insert = include_gxy + include_cs + include_drug
            insert2 = ''
    elif len(module_li)==1 and module_li[0]=='DRUG':
            insert = include_drug
            insert2 = ''
    else:
        if 'GXY' in module_li and 'CS' in module_li:
            if gxy=='1' and cs =='1':
                insert = include_gxy + include_cs + mingcishuoming
                insert2 = include_sanger
                if ordertag ==1: insert = include_gxy + include_cs + mingcishuoming + include_sanger
            elif gxy=='1' and cs =='0':
                insert = include_gxy + mingcishuoming + include_cs
                insert2 = include_sanger
                if ordertag ==1: insert = insert + include_sanger
            elif gxy=='0' and cs =='1':
                insert = include_gxy + include_cs + mingcishuoming
                insert2 = include_sanger
                if ordertag ==1: insert = insert + include_sanger
            elif gxy=='0' and cs =='0':
                insert = include_gxy + include_cs
                insert2 = ''
        elif 'GXY' not in module_li and 'CS' in module_li:
            if cs =='1':
                insert =  include_cs + mingcishuoming
                insert2 = include_sanger
                if ordertag ==1: insert = insert + include_sanger
            elif cs =='0':
                insert =  include_cs
                insert2 = ''
        elif 'CS' not in module_li and 'GXY' in module_li:
            if gxy=='1':
                insert = include_gxy + mingcishuoming
                insert2 = include_sanger
                if ordertag ==1: insert = insert + include_sanger
            elif gxy=='0':
                insert = include_gxy
                insert2 = ''
        if 'DRUG' in module_li:insert += include_drug
    CS_Dic['insert']=insert
    CS_Dic['insert2']=insert2
    return CS_Dic

def get_readfile(hospital):
    zfile='/PUBLIC/pipline/script/Report/XKFW/V3_180824/database/config.color.txt'
    covfile=''
    texfile=''
    fileDic={}
    if os.path.isfile(zfile) == False:
        print('%s is not a file' % zfile)
        sys.exit(1)
    try:
        with open(zfile, 'rU') as file_handle:
            for line in file_handle:
                line_list=[x.strip() for x in line.strip().split('\t')]
                fileDic[line_list[0]]=[line_list[1],line_list[2]]
            if '齐鲁' in hospital:
                covfile=fileDic['齐鲁'][0]
                texfile=fileDic['齐鲁'][1]
            else:
                covfile=fileDic['心康'][0]
                texfile=fileDic['心康'][1]
               
    except Exception, e:
        traceback.print_exc(e)
    return covfile,texfile

def main():
    argInit()
    if os.path.isfile(cs_file) == False:
        sys.stderr.write('%s is not a file' %cs_file)
        sys.exit(1)
    if os.path.isfile(gxy_file) == False:
        sys.stderr.write('%s is not a file' %gxy_file)
        sys.exit(1)
    cs_handle = open(cs_file,'r')
    cs_yinyang = cs_handle.readlines()[0].strip()
    cs_handle.close()
    gxy_handle = open(gxy_file,'r')
    gxy_yinyang = gxy_handle.readlines()[0].strip()
    gxy_handle.close()
    cover_Dic = defaultdict(list)
    sampleInf_Dic,hospital = get_sampleInf(sampleid,cfgFile,feature_file)
    coverfile,texfile=get_readfile(hospital)
    cover_Dic = get_Title_Letter(module,Module,coverfile)
    Shengming_Dic = get_Shengming(module)
    QC_Inf_Dic = get_QC_Inf(QCfile,sampleid)
    Jiance_Dic = get_Jiance(module)
    insert_Dic = get_insert_module(module,sampleid,gxy_yinyang,cs_yinyang,hospital)
    cover_Dic.update(sampleInf_Dic)
    cover_Dic.update(Shengming_Dic)
    cover_Dic.update(QC_Inf_Dic)
    cover_Dic.update(Jiance_Dic)
    cover_Dic.update(insert_Dic)
    out_all=open(outDir+'/'+sampleid+'_Cardio_all.tex','w')
    if cs_yinyang or gxy_yinyang:#wph add 20200413
        cover_Dic['Part'] = str(int(cover_Dic['Part'])+1) 
    with open(texfile,'r') as file_handle:
        cover_temp = Template(file_handle.read())
        report_cover = cover_temp.safe_substitute(cover_Dic)
        print >> out_all,report_cover
        #out_fengmian.write(report_cover)
    out_all.close()


if __name__ == '__main__':
    main()
