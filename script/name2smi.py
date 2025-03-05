# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 17:40:29 2019

@author: Administrator
"""

import sys
from ChemScript16 import StructureData

if __name__ == '__main__':
    name = sys.argv[1]
    m = StructureData()
    name = name.replace("_"," ")
    
    m.ReadData(name)
    smiles = m.Smiles
    print(smiles)