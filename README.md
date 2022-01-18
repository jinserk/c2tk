# C2TK: Computationl Chemistry Toolkit

[![DOI](https://zenodo.org/badge/414423138.svg)](https://zenodo.org/badge/latestdoi/414423138)

`c2tk` is a dockerized environment to run `psi4` with `qcengine`, `gpaw` and `orca` with `ase`:
- `python 3.9`
- source build:`libvdwxc-0.4.0`
- latest packages from pypi:
  - `python-dotenv`, `loguru`
  - `ase`, `gpaw`, `rdkit-pypi`
  - `networkx`, `matplotlib`, `nglview`
- install support for `orca-5.0.2` (if you download the binary from the [official site](https://orcaforum.kofo.mpg.de/app.php/dlext/))

## INSTALLATION

1. clone repo

2. download orca 5.0.2\
**ORCA is free only for academic use**\
register [here](https://orcaforum.kofo.mpg.de/index.php) and go downloads.\
click `ORCA 5.0.2` post and choose `ORCA 5.0.2, Linux, x86-64, shared-version, .tar.xz Archive`.
```
cd c2tk
mkdir .downloads
cp <your-download-path>/orca_5_0_2_linux_x86-64_shared_openmpi411.tar.xz .downloads
```

3. build image
```
./c2tk bulid
```

## HOW TO USE

Basically you can use docker and docker compose commands. We provide a script `./c2tk` for your convenience.

### Script help
```
./c2tk help
```

### Jupyter notebook
run `c2tk` container
```
./c2tk up
```

within ~1 minutes Jupyter notebook server will be ready.\
you can connect to the server via http://localhost:8989 .\
The default password is `Passw0rd!`.

### Shell environment
run `c2tk` shell
```
./c2tk bash
```

The default working dir is `/c2tk`.\
You can try to run `nwchem_openmpi` or `/opt/orca501/orca`.

### Python environment
run `ipython` shell
```
./c2tk ipython
```

or run your python scripts
```
./c2tk python examples/test.py
```

### Docker container down
```
./c2tk down
```

# CITATION
```
@software{Baik_C2TK_Computational_Chemistry_2021,
title = {{C2TK: Computational Chemistry Toolkit}},
author = {Baik, Jinserk},
doi = {10.5281/zenodo.5574254},
url = {https://github.com/jinserk/c2tk.git},
version = {0.0.1},
month = {10},
year = {2021}
}
```

# LICENSE
[BSD 3-Clause License](https://github.com/jinserk/c2tk/blob/main/LICENSE)
