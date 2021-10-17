# matk
Material Simulation Toolkit

[![DOI](https://zenodo.org/badge/414423138.svg)](https://zenodo.org/badge/latestdoi/414423138)

`matk` is a dockerized environment to run `nwchem`, `gpaw` and `orca` with `ase`:
- base OS image: `registry.fedoraproject.org/fedora-minimal:latest`
- `python 3.9.7`
- latest packages from fedora:
  - `environment-modules`, `openmpi`, `python3-devel`
  - `openblas-devel`, `openblas-openmp`, `openblas-threads`
  - `libxc-devel`, `fftw-devel`, `fftw-openmpi-devel`, `libomp-devel`
  - `elpa-openmpi-devel`, `blacs-openmpi-devel`, `hdf5-openmpi-devel`
  - `nwchem`, `nwchem-openmpi`
- source build:`libvdwxc-0.4.0`
- latest packages from pypi: `ase`, `gpaw`, `rdkit-pypi`, `networkx`, `matplotlib`
- install support for `orca-5.0.1` (if you download the binary from the [official site](https://orcaforum.kofo.mpg.de/app.php/dlext/))

## INSTALLATION

1. clone repo

2. download orca 5.0.1
**ORCA is free only for academic use**
register [here](https://orcaforum.kofo.mpg.de/index.php) and go downloads.\
click `ORCA 5.0.1` post and choose `ORCA 5.0.1, Linux, x86-64, shared-version, .tar.xz Archive`.\
```
cd matk
mkdir .downloads
cp <your-download-path>/orca_5_0_1_linux_x86-64_shared_openmpi411.tar.xz .downloads
```

3. build image
```
./matk bulid
```

## HOW TO USE

### Jupyter notebook
run `matk` container
```
./matk up
```

within ~1 minutes Jupyter notebook server will be ready.\
you can connect to http://localhost:8989. The default password is `Passw0rd!`

### Shell environment
run `matk` shell
```
./matk bash
```

The default working dir is `/matk`. You can try to run `nwchem_openmpi` or `/opt/orca501/orca`.

### Python environment
run `ipython` shell
```
./matk ipython
```

or run your python scripts
```
./matk python examples/main.py
```

### Docker container down
```
./matk down
```

# CITATION
```
@software{Baik_Material_Simulation_Toolki_2021,
author = {Baik, Jinserk},
doi = {10.5281/zenodo.5574253},
month = {10},
title = {{Material Simulation Toolki}},
url = {https://github.com/jinserk/matk},
version = {0.0.1},
year = {2021}
}
```

# LICENSE
[BSD 3-Cause License](https://github.com/jinserk/matk/blob/main/LICENSE)
