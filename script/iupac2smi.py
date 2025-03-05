# usage: python iupac2smi.py -i 'path to name csv file' -o 'path to name smiles map npy file'
import os,argparse
import numpy as np

name2cdx_script = './name2smi.py'
chemdraw_py32 = 'C:/Python32/python.exe'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, default='', help='path to name csv file')
    parser.add_argument('-o', '--output', type=str, default='', help='path to name smiles map npy file')
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    with open(input_file, 'r', encoding='utf-8') as fr:
        lines = fr.readlines()
    names = [line.strip() for line in lines]
    name_smiles_map = {}
    for idx,name in enumerate(names):
        name_blk_ = name.split(' ')
        name_blk = []
        for name_ in name_blk_:
            if not 'solution' in name_:
                name_blk.append(name_)
        new_name = '_'.join(name_blk)
        _ = os.popen('%s %s %s'%(chemdraw_py32,name2cdx_script,new_name))
        output_blks = _.read().split('\n')
        smi = output_blks[-2]
        print(f'[{idx+1} / {len(names)}] SMILES: {smi}')
        _.close()
        name_smiles_map[name] = smi
    np.save(output_file, name_smiles_map)
if __name__ == '__main__':
    main()
