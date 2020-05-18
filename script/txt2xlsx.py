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

def txt2xlx(source_txt,sheet_name="位点库",target_excel=None):
    if not target_excel:
        target_excel = os.path.splitext(source_txt)[0] + ".xlsx"
    wb = Workbook()
    sheet = wb.create_sheet("位点库".decode('utf-8'),index=0)
    source_file = open(source_txt,"r")
    for line in source_file.readlines():
        row_values = line.strip("\n").split("\t")
        for i in range(0,len(row_values)):
            if row_values[i].isdigit():
                row_values[i] = int(row_values[i])
        sheet.append(row_values)
    wb.save(target_excel)

if __name__ == "__main__":
    usage = """
    Usage:
        txt2xlsx.py <file> [-t <target_file>] [-n <sheet_name>]
        txt2xlsx.py -h
    
    Options:
        -h,--help                    显示帮助
        <file>                       输入txt文件
        -t,--target <target_file>    目标文件

    Optional:
        -n,--name <sheet_name>       sheet名字
        
"""
    from docopt import docopt
    args = docopt(usage)
    #print args
    source_txt = args['<file>']
    target_excel = args['--target'] 
    sheet_name = args['--name']
    print "txt2xlx " ,source_txt
    txt2xlx(source_txt,sheet_name,target_excel)
