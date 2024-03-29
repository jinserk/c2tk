#!/usr/bin/env python
from pathlib import Path

from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.qchem import QChem
from ase.optimize import LBFGS, QuasiNewton
from ase.vibrations import Vibrations

from gpaw import PW, FermiDirac

from c2tk.conformer import get_atoms, pre_optimize
from c2tk.calculators.nwchem import NWChem
from c2tk.calculators.orca import ORCA
from c2tk.calculators.gpaw import GPAW


def test_qchem(atoms: Atoms) -> None:
    calc = QChem(
        label='qchem',
        method='B3LYP',
        basis='6-31+G*',
        np=1, nt=4,
    )
    atoms.calc = calc

    #opt = LBFGS(atoms)
    #opt.run(fmax=0.05)
    e1 = atoms.get_potential_energy()
    print(f'test molecule energy: {e1:5.2f} eV')


def test_nwchem(atoms:Atoms) -> None:
    calc = NWChem(label='nwchem',
        dft=dict(
          maxiter=2000,
          xc='B3LYP',
        ),
        basis='6-31+G*',
    )
    atoms.calc = calc
    e1 = atoms.get_potential_energy()
    print(f'test molecule energy: {e1:5.2f} eV')

    """
    opt = LBFGS(atoms)
    opt.run(fmax=0.05)
    e2 = atoms.get_potential_energy()
    print(f'test molecule energy: {e2:5.2f} eV')

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
    print(f'test molecule energy: {e1:5.2f} eV')
    """
    relax = QuasiNewton(atoms, logfile='qn.log')
    relax.run(fmax=0.05)
    e2 = atoms.get_potential_energy()
    print(f'test molecule energy: {e2:5.2f} eV')
    """
    calc.write('temp.gpw')


def test_orca(atoms: Atoms) -> None:
    calc = ORCA(
        label='temp',
        orcasimpleinput='tightscf B3LYP/G def2-SVP kdiis opt freq',
        orcablocks='%scf maxiter 200 end\n%pal nprocs 8 end',
    )
    atoms.calc = calc
    e1 = atoms.get_potential_energy()
    print(f'test molecule energy: {e1:5.2f} eV')
    """
    relax = QuasiNewton(atoms, logfile='qn.log')
    relax.run(fmax=0.05)
    e2 = atoms.get_potential_energy()
    print(f'test molecule energy: {e2:5.2f} eV')
    """
    vib = Vibrations(atoms)
    vib.run()
    vib.summary()


def main(smiles: str) -> None:
    atoms = get_atoms(smiles)

    a1 = atoms.copy()
    a1.center(vacuum=50.0)

    a1.calc = EMT()
    e1 = a1.get_potential_energy()
    print(f'test molecule energy: {e1:5.2f} eV')

    """
    a2 = atoms.copy()
    a2.center(vacuum=50.0)
    a2 = pre_optimize(a2)

    a2.calc = EMT()
    e2 = a2.get_potential_energy()
    print(f'test molecule energy: {e2:5.2f} eV')
    """

    test_orca(atoms)
    #test_gpaw(atoms)
    #test_nwchem(atoms)


if __name__ == "__main__":
    smiles = "c1ccc2c(c1)[nH]c1ccc(-n3c4ccccc4c4ccccc43)cc12"
    main(smiles)

