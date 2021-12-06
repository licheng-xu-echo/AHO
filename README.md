# AHO
This is a repository for paper *"Towards Data-driven Design of **A**symmetric **H**ydrogenation of **O**lefins: Database and Hierarchical Learning"*. Here, you can find scripts used in this study.
# Introduction
Asymmetric hydrogenation of olefins is one of the most powerful asymmetric transformations in molecular synthesis. Although several privileged catalyst scaffolds are available, the catalyst development for asymmetric hydrogenation is still a time- and resource-consuming process due to the lack of predictive catalyst design strategy. Targeting the data-driven design of asymmetric catalysis, we herein report the development of a standardized database that contains the detailed information of over 12000 literature asymmetric hydrogenations of olefins. This database provides a valuable platform for the machine learning applications in asymmetric catalysis. Based on this database, we developed a hierarchical learning approach to achieve predictive machine leaning model using only dozens of enantioselectivity data with the target olefin, which offers a useful solution for the few-shot learning problem in the early stage of catalysis screening.
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
You can find information about reaction of Asymmetric Hydrogenation over [there](http://asymcatml.net/). Dataset availability is provided after registration
# How to cite
If the database or hierarchical learning is used in our publication, please cite as: Xu, L. -C.; Zhang, S. -Q.; Li, X.; Tang, M. -J.; Xie, P. -P.; Hong, X. *Angew. Chem. Int. Ed.* **2021**, *60*, 22804.
# Contact with us
Email: hxchem@zju.edu.cn; licheng_xu@zju.edu.cn
