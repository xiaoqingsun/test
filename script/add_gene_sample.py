#/bin/env/python
# -*- coding: utf-8 -*-
import sys
import xlrd
import os
reload(sys)
sys.setdefaultencoding('utf-8')
def SampleDiseaGene(sample_file,disea_gene_file,sample_disea_gene_file):
    sample_genes_dict = {}
    path = os.path.split(sample_disea_gene_file)
    disea_dict = ReadDisea(disea_gene_file)
    sample_disea_gene = open(sample_disea_gene_file,"w")
    (sample_list,sample_dict) = ReadCfg(sample_file)
    for sample in sample_list:
        if sample not in sample_genes_dict:
            sample_genes_dict[sample] = []
        #panel = ["心源性猝死","单基因高血压","单基因高血压三项","单基因糖尿病"]
        panel = ["单基因高血压三项"]
        disea_list = sample_dict[sample].split("+")
        if "高胆固醇血症" in disea_list:
            disea_list.append("谷固醇血症")
        if "低钾血症" in disea_list:
            disea_list.append("低血钾")
        if "心律失常" in disea_list:
            disea_list = disea_list + open("/PUBLIC/home/wangpenghui/心律失常").read().split("\n")
        for jibing in disea_list:
            if jibing in panel:
                continue
            elif jibing in "马凡综合征":
                jibing = "马方综合征"
            elif jibing.decode("utf-8") not in disea_dict:
                print "没有'%s'这个疾病的关联基因信息，或这个疾病名字有误，请检查\n" % jibing
            else:
                print jibing
                genes = disea_dict[jibing.decode("utf-8")]
                outline = "\t".join([sample,jibing,genes]) + "\n"
                if jibing in "单基因高血压":
                    genes = disea_dict["高血压".decode("utf-8")]
                    outline += "\t".join([sample,"高血压",genes]) + "\n"
                sample_disea_gene.write(outline)
                sample_genes_dict[sample] += genes.split(", ")
    sample_disea_gene.close()
    for sample in sample_list:
        gene_list_file = "%s/LOG/%s.genelist" %(path[0],sample)
        fp = open(gene_list_file,"w")
        for gene in set(sample_genes_dict[sample]): #wph changed tuple to set (to reduce)
            fp.write(gene + "\n")
        fp.close()
    return sample_disea_gene_file

def ReadCfg(cfg_file):
    cfg = open(cfg_file,"r")
    sample_dict = {}
    sample_list = []
    for line in cfg.read().split("\n")[1:-1]:
        sample = line.split("\t")[0]
        jibing = line.split("\t")[1]
        if "高脂血症" in jibing:
            jibing=jibing.replace("高脂血症","高胆固醇血症+高甘油三酯血症")
        sample_dict[sample] = jibing
        sample_list.append(sample)
    return sample_list,sample_dict


def ReadDisea(disea_gene_file):
    data = xlrd.open_workbook(disea_gene_file)
    table = data.sheets()[0]
    disea_dict = {}
    for i in range(table.ncols):
        gene_list = []
        disease = table.col(i)[1].value
        for j in range(2,len(table.col(i))):
            if "#" in table.col(i)[j].value:
                break
            else:
                gene_list.append(table.col(i)[j].value.strip())
        disea_dict[disease] = ", ".join(gene_list)
    return disea_dict

usage = """
    Python add_gene_sample.py report/01products/cfg.txt report/sample_disea_gene.txt
    """
if __name__ == "__main__":
    disea_gene_file = "/PUBLIC/pipline/database/Disease_WES/疾病-基因列表.xlsx"
    if len(sys.argv) ==1:
        print usage
    else:
        sample_file = sys.argv[1]
        sample_disea_gene_file = sys.argv[2]
        SampleDiseaGene(sample_file, disea_gene_file, sample_disea_gene_file)
