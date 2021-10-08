from ase.calculators.qchem import QChem
from ase.optimize import LBFGS

from matk.conformer import get_conformer
from matk.nwchem import NWChemWrapper, pre_optimize


def main(smiles: str) -> None:
    atoms = get_conformer(smiles)
    atoms.name = "test"

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


if __name__ == "__main__":
    smiles = "c1ccc2c(c1)[nH]c1ccc(-n3c4ccccc4c4ccccc43)cc12"
    main(smiles)
