import io

from ase import Atoms
from ase.io import xyz

from rdkit import Chem


class Molecule(Atoms):
    name = "anonymous"

    def to_xyz(self):
        with io.StringIO() as f:
            xyz.write_xyz(f, [self])
            return f.getvalue()
