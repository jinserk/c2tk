#!/usr/bin/env python

from ase import Atoms
from ase.calculators.qchem import QChem
from ase.calculators.nwchem import NWChem
from ase.calculators.orca import ORCA
from ase.optimize import LBFGS, QuasiNewton

from gpaw import GPAW, PW, FermiDirac

from matk.conformer import get_atoms, pre_optimize
#from matk.nwchem import NWChemWrapper, pre_optimize


def test_qchem(atoms: Atoms) -> None:
    calc = QChem(
        label='qchem',
        method='B3LYP',
        basis='6-31+G*',
        np=1, nt=4
    )
    atoms.calc = calc

    opt = LBFGS(atoms)
    opt.run(fmax=0.05)


def test_nwchem(atoms:Atoms) -> None:
    calc = NWChem(label='nwchem',
        dft=dict(
          maxiter=2000,
          xc='B3LYP',
        ),
        basis='6-31+G*'
    )
    atoms.calc = calc

    opt = LBFGS(atoms)
    opt.run(fmax=0.05)

    """
    obj = NWChemWrapper(nproc=1, mem=8000)
    calc_params = {
        'basis': '6-31+G*',
        'func': 'B3LYP',
        'target': 1,
    }
    e, f, p = obj.geom_opt(atoms, label="test", calc_params=calc_params)
    print(atoms.positions)
    """

def test_gpaw(atoms: Atoms) -> None:
    calc = GPAW(
        mode=PW(),
        xc='PBE',
        occupations=FermiDirac(0.0, fixmagmom=True),
        txt='temp.gpo',
    )
    atoms.calc = calc
    e1 = atoms.get_potential_energy()

    relax = QuasiNewton(atoms, logfile='qn.log')
    relax.run(fmax=0.05)
    e2 = atoms.get_potential_energy()

    print(f'hydrogen molecule energy: {e1:5.2f} eV')
    print(f'hydrogen molecule energy: {e2:5.2f} eV')
    calc.write('temp.gpw')


def test_orca(atoms: Atoms) -> None:
    calc = ORCA(
        label='temp',
        maxiter=2000, charge=0, mult=1, task='gradient',
        orcasimpleinput='B3LYP def2-SVP',
        orcablocks='%scf Convergence verytight\n maxiter 300 end\n %pal nprocs 1 end'
    )
    atoms.calc = calc
    e1 = atoms.get_potential_energy()

    relax = QuasiNewton(atoms, logfile='qn.log')
    relax.run(fmax=0.05)
    e2 = atoms.get_potential_energy()

    print(f'hydrogen molecule energy: {e1:5.2f} eV')
    print(f'hydrogen molecule energy: {e2:5.2f} eV')
    calc.write('temp.gpw')


def main(smiles: str) -> None:
    atoms = get_atoms(smiles)

    atoms.center(vacuum=5.0)
    atoms = pre_optimize(atoms)
    #test_orca(atoms)
    #test_gpaw(atoms)
    test_nwchem(atoms)


if __name__ == "__main__":
    smiles = "c1ccc2c(c1)[nH]c1ccc(-n3c4ccccc4c4ccccc43)cc12"
    main(smiles)
