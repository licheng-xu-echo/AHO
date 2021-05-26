# Some useful scripts

## 1. Scripts to convert .cdx files into SMILES and convert chemical compound name into .cdx
In this study, convert *.cdx* files into the correspond SMILES and in turn are one of the most time consuming work. So, we wrote the scripts to automatically interconvert the *.cdx* and SMILES.
In order to run *cdx2smi.py* and *iupac2cdx.py*, several third-party python packages and software are required.
```
ChemDraw
ChemScript16
comtypes>=1.1.9
pandas>=1.2.3
rdkit>=2019.09.3
```
The scripts mentioned above are executed under Windows 10. We use ChemDraw for easy to get started as it is widely used and other chemical structure drawing tools (e.g. ChemDoodle Web) should also work with neccessary modifications of these scripts.
### Usage
The function of *cdx2smi.py* is convert *.cdx* files into a *.csv* file that contains correspond SMILES.
```
python cdx2smil.py
please input the path which ".cdx" saved in (eg. C:/tmp/): **fill the right path here, then press Enter button**
please input one existed folder to save smiles file (eg. C:/tmp/smiles/): **fill the right path here, then press Enter button**
```
Here, we provide several *.cdx* files in *cdx_files* folder for test.  
  
The function of *iupac2cdx.py* is convert name of chemical compounds into correspond *.cdx* file and generate *.csv* file that contains SMILES information.
```
python iupac2cdx.py
please input the IUPAC name file: **fill the path of 'IUPAC_name.csv'**
please input the path of python.exe included in ChemDraw (default C:/Python32/python.exe):
please input the folder to save generated cdx files: **file the right path here**
```
Here, we provide a *IUPAC_name.csv* file in *name2cdx* folder for test.
