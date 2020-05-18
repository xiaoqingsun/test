#! /bin/env/python
# coding=utf-8

import os
import sys
import glob

def latest_conf():
    all_config_file_list = glob.glob("/PUBLIC/pipline/database/siteFilterDB/product/gxy_public_file/gxy_all_cfg_addsex/*.xls")
    latest_confi_file = ""
    last_time = 0
    for file in all_config_file_list:
        t = os.path.getmtime(file)
        if t > last_time:
            last_time = t
            latest_confi_file = file
    print latest_confi_file
    return latest_confi_file

def static_sample(latest_confi_file):
    sample_product_dict = {}
    sample_sex_dict = {}
    fp = open(latest_confi_file,"r")
    for line in fp.readlines()[1:]:
        items = line.strip("\n").split("\t")
        sample = items[0]
        sex = items[5]
        sample_sex_dict[sample] = sex
        pd_list = []
        if "三项" in items[1]:
            pd_list.append("GXY+CS")
        elif "高血压" in items[1] and "猝死" in items[1]:
            pd_list.append("GXY+CS")

        elif "高血压" in items[1] and "糖尿病" in items[1]:
            pd_list.append("GXY+DM")
        elif "猝死" in items[1] and "糖尿病" in items[1]:
            pd_list.append("CS+DM")
        elif "糖尿病" in items[1] or "唐尿病" in items[1]:
            pd_list.append("DM")
        elif "高血压" in items[1]:
            pd_list.append("GXY")
        elif "猝死" in items[1] or "心猝" in items[1]:
            pd_list.append("CS")
        elif "单基因疾病" in items[1]:
            pd_list.append("QR")
        elif "药物敏感" in items[1] or "用药指导" in items[1]:
            pd_list.append("Drug")
        elif "肿瘤检测" in items[1]:
            pd_list.append("Tumor")
        else:
            pd_list.append("WES")
    
        sample_product_dict[sample] = "+".join(pd_list)
    return sample_product_dict,sample_sex_dict
        #print sample,"\t",items[1],"\t",sample_product_dict[sample]



def split_sample_byproduct(sample_product_dict,samples):
    if samples == "-":
        return samples
    elif "." in samples:
        sample_list = samples.split(",")
        for i in range(len(sample_list)):
            if "." in sample_list[i]:
                sample_list[i] = sample_list[i].split(".")[0]
            elif "_" in sample_list[i]:
                sample_list[i] = sample_list[i].replace("_","")
        samples = ",".join(tuple(sample_list))
    sample_list = samples.split(",")
    result = ""
    product_list = ["GXY","CS","DM","WES","Drug","QR","Tumor"]
    product_dict = {}
    product_dict["GXY"] = []
    product_dict["CS"] = []
    product_dict["DM"] = []
    product_dict["WES"] = []
    product_dict["Drug"] = []
    product_dict["QR"] = []
    product_dict["Tumor"] = []
    for sample in sample_list:
        if "." in sample:
            sample = sample.split(".")[0]
        else:
            sample = sample.replace("_","")
        for product in sample_product_dict[sample].split("+"):
            product_dict[product].append(sample)
    for product in product_list:
        if len(product_dict[product]) == 0:
            del product_dict[product]
    for product in product_list:
        if product in product_dict:
            result += product + ":" + ",".join(tuple(product_dict[product])) + "|"
    return result[:-1]

def add_sample_product(sample_product_dict,sample_sex_dict,product,path):
    right_file = ""
    if product == "GXY":
        right_file = os.path.join(path,"GXY_sample_rigth_info.txt")
        wrong_file = os.path.join(path,"GXY_sample_wrong_info.txt")
        case_head = open(right_file,"r").readlines()[0].strip("\n").split("\t").index("case_ditaction_sample_ID")
        print case_head
    else:
        right_file = os.path.join(path,product+".right.txt")
        wrong_file = os.path.join(path,product+".err.txt")
        case_head = 92
    if "HET/HOM" in open(right_file,"r").readlines()[0]:
	het_col = open(right_file,"r").readlines()[0].strip("\n").split("\t").index("HET/HOM") 
    else:
        het_col = open(right_file,"r").readlines()[0].strip("\n").split("\t").index("Het/Hom")
    print "het_col:",het_col
    cmd = "cp %s %s.bak" %(right_file,right_file)
    os.system(cmd)
    cmd = "cp %s %s.bak" %(wrong_file,wrong_file)
    os.system(cmd)
    frigr = open(right_file + ".bak","r")
    frigw = open(right_file,"w")
    lines = frigr.readlines()
    
    headline = lines[0]
    frigw.write(headline)
    for line in lines[1:]:
        if "\t" not in line:
            frigw.write(line)
            continue
        items = line.strip("\n").split("\t")
        sample_id = items[1]
        chrom = items[2]
        if "X" in chrom and "男" in sample_sex_dict[sample_id]:
            items[het_col] = "hem"
        '''
        if product == "GXY" and items[8] == "NR3C2":
            split_list=[79,80,83,84,85,88,89,90,91]
            samples = items[case_head]
            if items[case_head] !="-":
                items[case_head] = split_sample_byproduct(sample_product_dict,samples)
            item1 = list(items)
            item2 = list(items)
            for i in split_list:
                if "|" in items[i]:
                    item1[i] = items[i].split("|")[0]
                    item2[i] = items[i].split("|")[1]
                    print item1[i],item2[i]
                else:
                    item1[i] = items[i]
                    item2[i] = items[i]
            outline = "\t".join(item1) + "\n" + "\t".join(item2) + "\n"
            frigw.write(outline)
            continue
        '''
        if len(items) < case_head:
            outline = "\t".join(items) + "\n"
            frigw.write(outline)
            continue
        samples = items[case_head]
        if samples == "-" or ":" in samples:
            outline = "\t".join(items) + "\n"
            frigw.write(outline)
            continue
        items[case_head] = split_sample_byproduct(sample_product_dict,samples)
        outline = "\t".join(items) + "\n"
        frigw.write(outline)
    frigr.close()
    frigw.close()

    fwronr = open(wrong_file + ".bak","r")
    fwronw = open(wrong_file,"w")
    for line in fwronr.readlines():
        items = line.strip("\n").split("\t")
        
        if len(items) < het_col or "case_ditaction_sample_ID" in line or "case_ID" in line:
            fwronw.write(line)
            continue
        sample_id = items[1]
        chrom = items[2]
        if "X" in chrom and "男" in sample_sex_dict[sample_id]:
            items[het_col] = "hem"
        if len(items) < case_head:
            outline = "\t".join(items) + "\n"
            fwronw.write(outline)
            continue
        samples = items[case_head]
        
        if samples == "-" or ":" in samples:
            outline = "\t".join(items) + "\n"
            fwronw.write(outline)
            continue
        items[case_head] = split_sample_byproduct(sample_product_dict,samples)
        outline = "\t".join(items) + "\n"
        fwronw.write(outline)
    fwronr.close()
    fwronw.close()
    

def main(args):
    product = args["--type"]
    path = args["--path"]
    latest_confi_file = latest_conf()
    sample_product_dict,sample_sex_dict = static_sample(latest_confi_file)
    add_sample_product(sample_product_dict,sample_sex_dict,product,path)
    


if __name__ == "__main__":
    usage = """
    Usage:
            add_sample_panel.py -t <product_type> -p <file_path>
            add_sample_panel.py -h
    
    Options:
            -h,--help                      显示帮助
            -t,--type <product_type>       分析产品类型(GXY,CS,DM,IDT)
            -p,--path <file_path>          right文件路径

    """
    from docopt import docopt
    args = docopt(usage)
    main(args)
