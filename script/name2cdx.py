# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 17:40:29 2019

@author: Administrator
"""

import sys
from ChemScript16 import StructureData

if __name__ == '__main__':
    name = sys.argv[1]
    dest_dir = sys.argv[2]
    m = StructureData()
    name = name.replace("_"," ")
    
    m.ReadData(name)
    smiles = m.Smiles
    m.ReadData(smiles)
    m.WriteFile('%s.cdx'%(dest_dir+'/'+name))
    print(smiles)