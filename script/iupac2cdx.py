# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 15:39:42 2021

@author: Administrator
"""


import pandas as pd
import os

## usage "python iupac2cdx.py"


if __name__ == '__main__':
    name2cdx_script = './name2cdx.py'
    iupac_file_name = input('please input the IUPAC name file: ')
    
    default_py32 = 'C:/Python32/python.exe'
    py32 = input('please input the path of python.exe included in ChemDraw (default C:/Python32/python.exe): ')
    if py32 == '':
        py32 = default_py32
    cdx_dir = input('please input the folder to save generated cdx files: ')
    
    df = pd.read_csv(iupac_file_name,index_col=0)
    
    name = df['IUPAC name'].to_list()
    name = [item.replace(" ","_") for item in name]
    iupac_smi_dict = {'IUPAC name':[],'SMILES':[]}
    for item in name:
        iupac_smi_dict['IUPAC name'].append(item)
        _ = os.popen('%s %s %s %s'%(py32,name2cdx_script,item,cdx_dir))
        smi = _.read().split('\n')[-2]
        _.close()
        iupac_smi_dict['SMILES'].append(smi)
    df = pd.DataFrame.from_dict(iupac_smi_dict)
    df.to_csv(cdx_dir+'/iupac_smiles.csv')
