import io

import ase
from ase.io import extxyz
from ase.calculators.emt import EMT

from rdkit import Chem

from .xyz2mol import xyz2mol

class Atoms(ase.Atoms):
    name = "anonymous"

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
