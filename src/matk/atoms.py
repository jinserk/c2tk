import io
import random

import ase
from ase.io import extxyz

from rdkit import Chem
from rdkit.Chem.AllChem import (
    EmbedMultipleConfs, GetBestRMS,
    MMFFGetMoleculeForceField, MMFFGetMoleculeProperties,
    UFFGetMoleculeForceField,
)

from .conformer import get_atoms, pre_optimize
from .xyz2mol import xyz2mol


class Atoms(ase.Atoms):

    def add_mol(self, smiles: str,
                center: list[float] = None,
                rotate: list[float] = None,
                pre_opt: bool = False):
        atoms = get_atoms(smiles)
        if pre_opt:
            pre_optimize(atoms)
        self.extend(atoms)

    def to_rdmol(self):
        atoms = self.get_atomic_numbers().tolist()
        xyz_pos = [(p[0], p[1], p[2]) for p in self.positions]
        charge = int(sum(self.get_initial_charges()))
        mols = xyz2mol(atoms, xyz_pos, charge=charge,
                       use_graph=True,
                       allow_charged_fragments=True,
                       embed_chiral=True,
                       use_huckel=True)
        return mols

    def to_xyz(self):
        with io.StringIO() as f:
            extxyz.write_extxyz(f, [self])
            return f.getvalue()


