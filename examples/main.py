#from ase.calculators.qchem import QChem
from gpaw import GPAW, PW, FermiDirac
from ase.optimize import QuasiNewton

from matk.conformer import get_atoms, pre_optimize
#from matk.nwchem import NWChemWrapper, pre_optimize


'''
def main(smiles: str) -> None:
    atoms = get_conformer(smiles)
    atoms.name = "test"
    print(atoms.to_rdmol())

    atoms = pre_optimize(atoms)

    """
    calc = QChem(label='qchem',
                 method='B3LYP',
                 basis='6-31+G*',
                 np=1, nt=4)
    atoms.calc = calc

    opt = LBFGS(atoms)
    opt.run(fmax=0.05)

    calc = NWChem(label='nwchem',
                  dft=dict(
                      maxiter=2000,
                      xc='B3LYP',
                  ),
                  basis='6-31+G*')
    calc.command = f"mpirun -np {nproc} nwchem PREFIX.nwi > PREFIX.nwo"
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
'''

def main(smiles: str) -> None:
    atoms = get_atoms(smiles)

    atoms.center(vacuum=5.0)
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

if __name__ == "__main__":
    smiles = "c1ccc2c(c1)[nH]c1ccc(-n3c4ccccc4c4ccccc43)cc12"
    main(smiles)
