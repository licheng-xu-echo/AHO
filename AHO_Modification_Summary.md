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
