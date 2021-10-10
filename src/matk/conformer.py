import random

from rdkit import Chem
from rdkit.Chem.AllChem import (
    EmbedMultipleConfs, GetBestRMS,
    MMFFGetMoleculeForceField, MMFFGetMoleculeProperties,
    UFFGetMoleculeForceField,
)

from .atoms import Atoms
from .nwchem import NWChemWrapper


def make_atoms(mol, attempts=50):
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

    aa = []
    pp = []
    atoms = mol.GetAtoms()
    for atom_num in range(0, conf.GetNumAtoms()):
        pos = conf.GetAtomPosition(atom_num)
        a = atoms[atom_num].GetSymbol()
        aa.append(a)
        pp.append((pos.x, pos.y, pos.z))

    return Atoms(aa, positions=pp)


def get_conformer(smiles: str) -> Atoms:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    atoms = make_atoms(mol)
    return atoms
