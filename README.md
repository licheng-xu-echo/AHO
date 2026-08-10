# AHO
This is a repository for paper [*"Towards Data-driven Design of **A**symmetric **H**ydrogenation of **O**lefins: Database and Hierarchical Learning"*](https://onlinelibrary.wiley.com/doi/10.1002/anie.202106880). Here, you can find scripts used in this study.

👉️**The full AHO dataset is available at [AHO-Dataset.zip](https://github.com/licheng-xu-echo/AHO/blob/main/AHO-Dataset.zip).** 👈️

The online database server can be accessed at http://asymcatml.net/
# Introduction
Asymmetric hydrogenation of olefins is one of the most powerful asymmetric transformations in molecular synthesis. Although several privileged catalyst scaffolds are available, the catalyst development for asymmetric hydrogenation is still a time- and resource-consuming process due to the lack of predictive catalyst design strategy. Targeting the data-driven design of asymmetric catalysis, we herein report the development of a standardized database that contains the detailed information of over 12000 literature asymmetric hydrogenations of olefins. This database provides a valuable platform for the machine learning applications in asymmetric catalysis. Based on this database, we developed a hierarchical learning approach to achieve predictive machine leaning model using only dozens of enantioselectivity data with the target olefin, which offers a useful solution for the few-shot learning problem in the early stage of catalysis screening.
![AHO_TOC](meta/anie202106880-toc-0001-m.jpg)
# Dependence
In order to run Jupyter Notebook for machine learning application demonstration, several third-party python packages are required.
```
python>=3.8.5
numpy>=1.19.2
pandas>=1.2.0
ase>=3.21.0
dscribe>=1.0.0
rdkit>=2019.09.3
openbabel>=3.1.0
scikit-learn>=0.23.2
mordred>=1.2.0
matplotlib>=3.3.2
```
We suggest using [Anaconda](https://www.anaconda.com/) to install the python 3.8.5 or higher version, as conda and pip together make the installation of these dependences much easier. All test are executed under Ubuntu 18.04, as the [*dscribe*](https://singroup.github.io/dscribe/latest/install.html) package currently only support Unix-based systems.
# Installation of dependence
We suggest using [Anaconda](https://www.anaconda.com/) to prepare dependence as many packages are built-in Anaconda base environment. For those packages not built-in, you may input following commands to install them and follow the installation instructions.
```
conda install ase
pip install dscribe
conda install rdkit -c rdkit
conda install openbabel -c conda-forge
conda install -c rdkit -c mordred-descriptor mordred
```
# Usage
Here we provide [several tutorials](https://github.com/licheng-xu-echo/AHO/tree/main/examples) in Jupyter Notebook format to demonstrate how to generate descriptors with provided reaction data, train machine learning model and use *hierarchical learning* approach to handle few-shot learning problem.
# Dataset availability
You can find information about reaction of Asymmetric Hydrogenation over [there](http://asymcatml.net/).
# How to cite
If the database or hierarchical learning is used in your publication, please cite as: Xu, L. -C.; Zhang, S. -Q.; Li, X.; Tang, M. -J.; Xie, P. -P.; Hong, X. *Angew. Chem. Int. Ed.* **2021**, *60*, 22804.
# Contact with us
Email: hxchem@zju.edu.cn; licheng_xu@zju.edu.cn

---

# 📝 AHO Modification Change Log

**Project Name**: Asymmetric Hydrogenation of Olefins (AHO)  
**Goal**: Fix Windows compatibility, resolve expired external dataset URLs, and improve code stability under modern Matplotlib and RDKit environments.  
**Test Environment**: Windows 11/10 | Python 3.10 | RDKit | Scikit-Learn | Matplotlib  

---

## 📂 Summary of Modified Files

### 1. Core Script & Utility Modifications

- **`examples/mlutils.py`**  
  Fixed invalid LaTeX math syntax in `drawregfig` by replacing `\itG` with standard `G` (ΔΔG). Added `try-except` exception guard around `plt.tight_layout()`.

---

### 2. Tutorial Notebook Fixes

- **`examples/01-generate_descriptors.ipynb`**  
  Replaced Linux-specific bash commands (`! wget`, `! unzip`) with cross-platform Python `zipfile` extraction logic.

- **`examples/02-train_simple_ML_model.ipynb`**  
  Replaced `%matplotlib notebook` with `%matplotlib inline`. Added fallback logic to extract local `datasets_for_fig_4.zip` and generate dataset matrices (resolving expired `asymcatml.net` URLs).

- **`examples/03-hierarchical_learning.ipynb`**  
  Wrapped `from rdkit.Chem import Draw` in `try-except` to prevent Windows GUI DLL load failures. Replaced `! cat /proc/cpuinfo` with `os.cpu_count()`. Updated `display.set_matplotlib_formats('svg')` to standard `set_matplotlib_formats('svg')`.

- **`examples/04-generate_correct_chiral_biaryl.ipynb`**  
  Added `os.makedirs('./biaryl', exist_ok=True)` to ensure `mol_1.sdf` exports smoothly.

- **`examples/05-generate_correct_molecule_contain_ferrocene.ipynb`**  
  Added `try-except` guard for `Draw` import and display. Included missing ferrocene template `examples/Ferr/Ferr_std_reverse.sdf`.

- **`examples/06-show_olefins_approaching_PCA.ipynb`**  
  Replaced `%matplotlib notebook` with `%matplotlib inline`. Added `Draw` import guard.

---

### 3. Added Template & Preprocessed Datasets

- **`examples/Ferr/Ferr_std_reverse.sdf`**  
  Generated reverse ferrocene 3D template file by flipping X coordinates of `Ferr_std.sdf`.

- **`examples/data/screening_desc_ensemble.npz`**  
  Preprocessed descriptor ensemble dataset generated from local `related_dataset_a.csv`.

- **`examples/data/target.npz`**  
  Preprocessed target dataset.

- **`examples/data/hierarchical_learning_set.npz`**  
  Preprocessed hierarchical learning dataset.

---

## 🛠️ Key Improvement Areas

### 1. Cross-Platform Compatibility (Windows vs Linux)
- **Problem**: Original code used Linux-specific bash commands (`! wget`, `! unzip`, `! mv`, `! cat /proc/cpuinfo`), causing syntax errors or missing command failures on Windows.
- **Solution**: Refactored logic using Python standard library modules (`os`, `zipfile`, `sys`) to ensure out-of-the-box execution across Windows, macOS, and Linux.

### 2. Expired External Data Source Fallback
- **Problem**: Original code relied on external laboratory server `http://asymcatml.net/download/...`, which is now unreachable (404 Error).
- **Solution**: Added automated extraction and RDKit descriptor fallback using local `datasets_for_fig_4.zip`, eliminating external server dependency.

### 3. Modern Matplotlib & RDKit Syntax Compatibility
- **Problem**: 
  - Matplotlib math parser throws `ValueError: Unknown symbol: \itG` for legacy `\itG` LaTeX syntax.
  - `%matplotlib notebook` triggers `IPython is not defined` JavaScript error in modern JupyterLab.
  - RDKit `Draw` import fails on Windows environments missing GUI/Cairo C++ DLLs.
- **Solution**: Switched to `%matplotlib inline`, updated LaTeX `\itG` to `G`, and added `try-except` guards around `Draw` imports.
