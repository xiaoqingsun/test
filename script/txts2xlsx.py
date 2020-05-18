#!/usr/bin/env python2
#coding:utf-8
import os
import sys
import datetime
import time
import openpyxl
from openpyxl import Workbook
reload(sys)
sys.setdefaultencoding("utf-8")

def txt2xlx(source_txt,wb,sheet_name="位点库",sheet_index=0):
    sheet = wb.create_sheet(sheet_name,index=sheet_index)
    source_file = open(source_txt,"r")
    for line in source_file.readlines():
        row_values = line.strip("\n").split("\t")
        for i in range(0,len(row_values)):
            if row_values[i].isdigit():
                row_values[i] = int(row_values[i])
        sheet.append(row_values)

def MutiFiles(source_txt,sheet_name="位点库",target_excel=None):
    if not target_excel:
        target_excel = os.path.splitext(source_txt)[0] + ".xlsx"
    wb = Workbook()
    txts=source_txt.split(",")
    names=[]
    if sheet_name:
        names=sheet_name.replace("，",",").split(",")
    for i in range(len(txts)):
        if not txts[i]:
            continue
        if len(names) <i or len(names)==0:
            sheet="sheet%s" %i
        else:
            sheet=names[i]
        print sheet
        txt2xlx(txts[i],wb,unicode(sheet),i)
    wb.save(target_excel)


if __name__ == "__main__":
    usage = """
    Usage:
        txts2xlsx.py <file> [-t <target_file>] [-n <sheet_name>]
        txts2xlsx.py -h
    
    Options:
        -h,--help                    显示帮助
        <file>                       输入txt文件，多个文件用“,”分隔
        -t,--target <target_file>    目标文件

    Optional:
        -n,--name <sheet_name>       sheet名字,多个sheet用“,”分隔
        
"""
    from docopt import docopt
    args = docopt(usage)
    source_txt = args['<file>']
    target_excel = args['--target'] 
    sheet_name = args['--name']
    print "txts2xlx " ,source_txt
    MutiFiles(source_txt,sheet_name,target_excel)
