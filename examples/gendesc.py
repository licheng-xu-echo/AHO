# -*- coding: utf-8 -*-
"""


@author: Li-Cheng Xu
"""
import numpy as np
import glob
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from mordred import Calculator, descriptors
from dscribe.descriptors import ACSF,SOAP,LMBTR,MBTR
from ase import Atoms as ASE_Atoms
period_table = Chem.GetPeriodicTable()

def getmorganfp(mol,radius=2,nBits=2048,useChirality=True):
    '''
    
    Parameters
    ----------
    mol : mol
        RDKit mol object.

    Returns
    -------
    mf_desc_map : ndarray
        ndarray of molecular fingerprint descriptors.

    '''
    fp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,radius=radius,nBits=nBits,useChirality=useChirality)
    return np.array(list(map(eval,list(fp.ToBitString()))))


def Mol2Atoms(mol):
    positions = mol.GetConformer().GetPositions()
    atom_types = [period_table.GetElementSymbol(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
    atoms = ASE_Atoms(symbols=atom_types,positions=positions)
    return atoms



def getusidx(desc):
    usidx = []
    desc_scale = desc.max(axis=0) - desc.min(axis=0)
    for i in range(len(desc_scale)):
        if np.isnan(desc_scale[i]) or desc_scale[i] == 0:
            continue
        usidx.append(i)
    return usidx
class generate2Ddesc():

    def __init__(self,df):
        
        '''
        
        Parameters
        ----------
        df : DataFrame(pandas)
            dataframe which needed process.

        Returns
        -------
        gen2d: obeject
        '''
        
        self.re_smi = df['Reactant SMILES'].to_numpy()
        self.pr_smi = df['Product SMILES'].to_numpy()
        self.sol_smi = df['Solvent SMILES'].to_numpy()
        self.cat_smi = df['Catalyst SMILES(RDKit)'].to_numpy()

        self.smi_set = np.unique(np.concatenate([self.re_smi,self.pr_smi,self.sol_smi,self.cat_smi]))
        self.mol_set = [Chem.MolFromSmiles(tmp_smi) for tmp_smi in self.smi_set]
    def getmorganfp(self,mol,radius=2,nBits=2048,useChirality=True):
        '''
        
        Parameters
        ----------
        mol : mol
            RDKit mol object.

        Returns
        -------
        mf_desc_map : ndarray
            ndarray of molecular fingerprint descriptors.

        '''
        fp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,radius=radius,nBits=nBits,useChirality=useChirality)
        return np.array(list(map(eval,list(fp.ToBitString()))))
    def calc_rdkit_desc(self):
        '''
        
        Parameters
        ----------
        

        Returns
        -------
        rdkit_desc_map : dict
            map of RDKit descriptors.

        '''
        descs = [desc_name[0] for desc_name in Descriptors._descList]
        desc_calc = MoleculeDescriptors.MolecularDescriptorCalculator(descs)
        rdkit_desc_map = {self.smi_set[i] : np.array(desc_calc.CalcDescriptors(self.mol_set[i])) for i in range(len(self.smi_set))}
        return rdkit_desc_map
    def calc_modred_desc(self):
        '''
        

        Returns
        -------
        modred_desc_map : dict
            map of modred descriptors.

        '''
        calc = Calculator(descriptors, ignore_3D=True)
        modred_df = calc.pandas(self.mol_set)
        modred_desc_map = {self.smi_set[i]: modred_df.iloc[i].to_numpy() for i in range(len(self.smi_set))}
        return modred_desc_map
    def calc_mf_desc(self):
        '''
        

        Returns
        -------
        mf_desc_map : dict
            map of molecular fingerprint descriptors.

        '''
        
        mf_desc_map = {self.smi_set[i]:self.getmorganfp(self.mol_set[i]) for i in range(len(self.smi_set))}
        return mf_desc_map

class generate3Ddesc():
    def __init__(self,sdf_file_dir):
        self.sdf_files = glob.glob(sdf_file_dir+'*.sdf')
        self.metal_types = ['Rh','Ir','Co','Ru','Ni','Pd']
        #xyz_files = glob.glob('./data/geoms_demo/*.xyz')
    def getkeyatom(self,tmp_mol):
        pos = tmp_mol.GetConformer().GetPositions()
        atom_weights = np.array([tmp_atom.GetMass() for tmp_atom in tmp_mol.GetAtoms()]).reshape(-1,1)
        atom_weights = np.concatenate([atom_weights,atom_weights,atom_weights],axis=1)
        weight_cent = np.sum(pos*atom_weights,axis=0)/atom_weights.sum()
        key_atom = np.argmin(np.sum((pos - weight_cent)**2,axis=1))
        return key_atom
    def Mol2Atoms(self,mol):
        positions = mol.GetConformer().GetPositions()
        atom_types = [period_table.GetElementSymbol(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
        atoms = ASE_Atoms(symbols=atom_types,positions=positions)
        return atoms
    def getkeyatomspecies(self):
        
        key_atoms = []
        atom_species = {'react':[],'prod':[],'sol':[],'cat':[]}
        
        for tmp_sdf_file in self.sdf_files:
            with open(tmp_sdf_file,'r') as fr:
                type_ = fr.readlines()[0].strip().split()[1]
            tmp_mol = Chem.MolFromMolFile(tmp_sdf_file,sanitize=False,removeHs=False)
            tmp_sym = [tmp_at.GetSymbol() for tmp_at in tmp_mol.GetAtoms()]
            try:
                atom_species[type_] += tmp_sym
            except:
                atom_species['cat'] += tmp_sym
            key_atoms.append(self.getkeyatom(tmp_mol))
        atom_species = {tmp_key:list(set(atom_species[tmp_key])) for tmp_key in atom_species}
        self.key_atoms = key_atoms
        self.atom_species = atom_species
        return key_atoms,atom_species
    
    def calc_acsf_desc(self):
        key_atoms,atom_species = self.getkeyatomspecies()
        rcut=6.0
        g2_params=[[1, 1], [1, 2], [1, 3]]
        g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1]]
        cat_calc = ACSF(species=atom_species['cat'],rcut=rcut,g2_params=g2_params,g4_params=g4_params)
        re_calc = ACSF(species=atom_species['react'],rcut=rcut,g2_params=g2_params,g4_params=g4_params)
        pr_calc = ACSF(species=atom_species['prod'],rcut=rcut,g2_params=g2_params,g4_params=g4_params)
        sol_calc = ACSF(species=atom_species['sol'],rcut=rcut,g2_params=g2_params,g4_params=g4_params)
        acsf_desc_map = {}
        for i in range(len(self.sdf_files)):
            tmp_sdf_file = self.sdf_files[i]
            with open(tmp_sdf_file,'r') as fr:
                title_line = fr.readlines()[0].strip().split()
                tmp_smi = title_line[0]
                type_ = title_line[1]
            tmp_mol = Chem.MolFromMolFile(tmp_sdf_file,sanitize=False,removeHs=False)
            tmp_atom = self.Mol2Atoms(tmp_mol)

            if type_ == 'react':
                tmp_desc = re_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'prod':
                tmp_desc = pr_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'sol':
                tmp_desc = sol_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            else:
                atom_syms = [tmp_atom.GetSymbol() for tmp_atom in tmp_mol.GetAtoms()]
                for tmp_m in self.metal_types:
                    try:
                        metal_key_atom_idx = atom_syms.index(tmp_m)
                        break
                    except:
                        continue
                tmp_desc = cat_calc.create(tmp_atom,positions=[metal_key_atom_idx,key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            acsf_desc_map[tmp_smi] = tmp_desc
        return acsf_desc_map
    
    def calc_soap_desc(self):
        key_atoms,atom_species = self.getkeyatomspecies()
        rcut = 6.0
        nmax = 4
        lmax = 3
        re_calc = SOAP(species=atom_species['react'],rcut=rcut,nmax=nmax,lmax=lmax)
        pr_calc = SOAP(species=atom_species['prod'],rcut=rcut,nmax=nmax,lmax=lmax)
        cat_calc = SOAP(species=atom_species['cat'],rcut=rcut,nmax=nmax,lmax=lmax)
        sol_calc = SOAP(species=atom_species['sol'],rcut=rcut,nmax=nmax,lmax=lmax)
        soap_desc_map = {}
        for i in range(len(self.sdf_files)):
            tmp_sdf_file = self.sdf_files[i]
            with open(tmp_sdf_file,'r') as fr:
                title_line = fr.readlines()[0].strip().split()
                tmp_smi = title_line[0]
                type_ = title_line[1]
            tmp_mol = Chem.MolFromMolFile(tmp_sdf_file,sanitize=False,removeHs=False)
            tmp_atom = self.Mol2Atoms(tmp_mol)

            if type_ == 'react':
                tmp_desc = re_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'prod':
                tmp_desc = pr_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'sol':
                tmp_desc = sol_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            else:
                atom_syms = [tmp_atom.GetSymbol() for tmp_atom in tmp_mol.GetAtoms()]
                for tmp_m in self.metal_types:
                    try:
                        metal_key_atom_idx = atom_syms.index(tmp_m)
                        break
                    except:
                        continue
                tmp_desc = cat_calc.create(tmp_atom,positions=[metal_key_atom_idx,key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            soap_desc_map[tmp_smi] = tmp_desc
        return soap_desc_map
    def calc_lmbtr_desc(self):
        key_atoms,atom_species = self.getkeyatomspecies()
        k2={
                "geometry": {"function": "inverse_distance"},
                "grid": {"min": 0, "max": 1, "n": 10, "sigma": 0.1},
                "weighting": {"function": "exponential", "scale": 0.5, "cutoff": 1e-3},
            }
        k3={
                "geometry": {"function": "cosine"},
                "grid": {"min": -1, "max": 1, "n": 10, "sigma": 0.1},
                "weighting": {"function": "exponential", "scale": 0.5, "cutoff": 1e-3},
            }
        re_calc = LMBTR(species=atom_species['react'],k2=k2,k3=k3,periodic=False)
        pr_calc = LMBTR(species=atom_species['prod'],k2=k2,k3=k3,periodic=False)
        cat_calc = LMBTR(species=atom_species['cat'],k2=k2,k3=k3,periodic=False)
        sol_calc = LMBTR(species=atom_species['sol'],k2=k2,k3=k3,periodic=False)
        lmbtr_desc_map = {}
        for i in range(len(self.sdf_files)):
            tmp_sdf_file = self.sdf_files[i]
            with open(tmp_sdf_file,'r') as fr:
                title_line = fr.readlines()[0].strip().split()
                tmp_smi = title_line[0]
                type_ = title_line[1]
            tmp_mol = Chem.MolFromMolFile(tmp_sdf_file,sanitize=False,removeHs=False)
            tmp_atom = self.Mol2Atoms(tmp_mol)

            if type_ == 'react':
                tmp_desc = re_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'prod':
                tmp_desc = pr_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'sol':
                tmp_desc = sol_calc.create(tmp_atom,positions=[key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            else:
                atom_syms = [tmp_atom.GetSymbol() for tmp_atom in tmp_mol.GetAtoms()]
                for tmp_m in self.metal_types:
                    try:
                        metal_key_atom_idx = atom_syms.index(tmp_m)
                        break
                    except:
                        continue
                tmp_desc = cat_calc.create(tmp_atom,positions=[metal_key_atom_idx,key_atoms[i]])
                tmp_desc = np.concatenate(tmp_desc)
            lmbtr_desc_map[tmp_smi] = tmp_desc
        return lmbtr_desc_map
    def calc_mbtr_desc(self):
        key_atoms,atom_species = self.getkeyatomspecies()
        k1={
                "geometry": {"function": "atomic_number"},
                "grid": {"min": 0, "max": 8, "n": 10, "sigma": 0.1},
            }
        k2={
                "geometry": {"function": "inverse_distance"},
                "grid": {"min": 0, "max": 1, "n": 10, "sigma": 0.1},
                "weighting": {"function": "exponential", "scale": 0.5, "cutoff": 1e-3},
            }
        k3={
                "geometry": {"function": "cosine"},
                "grid": {"min": -1, "max": 1, "n": 10, "sigma": 0.1},
                "weighting": {"function": "exponential", "scale": 0.5, "cutoff": 1e-3},
            }
        re_calc = MBTR(species=atom_species['react'],k1=k1,k2=k2,k3=k3,periodic=False,normalization="l2_each")
        pr_calc = MBTR(species=atom_species['prod'],k1=k1,k2=k2,k3=k3,periodic=False,normalization="l2_each")
        cat_calc = MBTR(species=atom_species['cat'],k1=k1,k2=k2,k3=k3,periodic=False,normalization="l2_each")
        sol_calc = MBTR(species=atom_species['sol'],k1=k1,k2=k2,k3=k3,periodic=False,normalization="l2_each")
        mbtr_desc_map = {}
        for i in range(len(self.sdf_files)):
            tmp_sdf_file = self.sdf_files[i]
            with open(tmp_sdf_file,'r') as fr:
                title_line = fr.readlines()[0].strip().split()
                tmp_smi = title_line[0]
                type_ = title_line[1]
            tmp_mol = Chem.MolFromMolFile(tmp_sdf_file,sanitize=False,removeHs=False)
            tmp_atom = self.Mol2Atoms(tmp_mol)
            
            if type_ == 'react':
                tmp_desc = re_calc.create(tmp_atom)
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'prod':
                tmp_desc = pr_calc.create(tmp_atom)
                tmp_desc = np.concatenate(tmp_desc)
            elif type_ == 'sol':
                tmp_desc = sol_calc.create(tmp_atom)
                tmp_desc = np.concatenate(tmp_desc)
            else:
                tmp_desc = cat_calc.create(tmp_atom)
                tmp_desc = np.concatenate(tmp_desc)
            mbtr_desc_map[tmp_smi] = tmp_desc
        return mbtr_desc_map
    