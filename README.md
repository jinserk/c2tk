# matk
Material Simulation Toolkit

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

## Installation

### clone repo

### download orca 5.0.1
register [here](https://orcaforum.kofo.mpg.de/index.php) and go downloads.\
click `ORCA 5.0.1` post and choose `ORCA 5.0.1, Linux, x86-64, shared-version, .tar.xz Archive`.\
make a folder named `.downloads/` under the `matk/`
put your downloaded xz file to `matk/.downloads/

### build image
```
./matk bulid
```

## How to use

### Jupyter notebook
run `matk` container
```
./matk up
```

after 1 or 2 min, you can connect to http://localhost:8989. The default password is `Passw0rd!`

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
