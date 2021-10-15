import random
from typing import Optional

import ase
from ase.calculators.emt import EMT
from ase.optimize import LBFGS

from rdkit import Chem
from rdkit.Chem.AllChem import (
    EmbedMultipleConfs, GetBestRMS,
    MMFFGetMoleculeForceField, MMFFGetMoleculeProperties,
    UFFGetMoleculeForceField,
)


def guess_conformer(mol: Chem.rdchem.Mol,
                    attempts: int = 50) -> Chem.rdchem.Conformer:
    max_attempts = attempts * 25
    EmbedMultipleConfs(mol,
                       numConfs=attempts,
                       maxAttempts=max_attempts,
                       useRandomCoords=True,
                       randomSeed=random.randint(1, 10000000))

    num_confs = len(mol.GetConformers())
    if not num_confs:
        raise ValueError(f"No conformer with {max_attempts} attemps")

    conf_energies = []
    for i in range(0, num_confs):
        try:
            props = MMFFGetMoleculeProperties(mol)
            potential = MMFFGetMoleculeForceField(mol, props, confId=i)
            if potential is None:
                potential = UFFGetMoleculeForceField(mol, confId=i)
            potential.Minimize()
            ff_energy = potential.CalcEnergy()
            conf_energies.append((i, ff_energy))
        except:
            continue

    min_id = sorted(conf_energies, key=lambda x: x[1])[0][0]
    conf = mol.GetConformer(id=min_id)
    return conf


def make_atoms(mol: Chem.rdchem.Mol) -> ase.Atoms:
    conf = guess_conformer(mol)
    mol_atoms = mol.GetAtoms()

    atoms = ase.Atoms()
    for an in range(0, conf.GetNumAtoms()):
        a = mol_atoms[an].GetSymbol()
        p = conf.GetAtomPosition(an)
        atoms.append(ase.Atom(a, [p.x, p.y, p.z]))
    return atoms


def get_atoms(smiles: str) -> ase.Atoms:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    atoms = make_atoms(mol)
    return atoms


def pre_optimize(atoms: ase.Atoms,
                 fmax: float = 0.05) -> None:
    _calc = atoms.calc
    atoms.calc = EMT()
    opt = LBFGS(atoms, logfile=f"temp.pre")
    opt.run(fmax=fmax)
    atoms.calc = _calc
    return atoms

