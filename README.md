DocKin
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/DocKin.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/DocKin)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/DocKin/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/DocKin/branch/master)

A Python library for docking small molecules into kinases.

The `DocKin` repository is used to play with different docking algorithms and to investigate their performance with focus on kinases.

Integrated docking tools:

- [OEDock](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) (license required)

### Install

1. Clone repository

`git clone https://github.com/volkamerlab/DocKin.git`

2. Create Conda environment

`cd devtools/conda-envs`  
`conda env create -f env.yml`  
`conda activate dockin`

### Examples

Jupyter notebooks are provided in the example directory providing a step-by-step guide on how to use implemented functions.

### Copyright

Copyright (c) 2020, Volkamer Lab

### License

`DocKin` is free software and is licensed under the MIT license.

### Authors

- David Schaller <david.schaller@charite.de>

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.2.
