from pathlib import Path
import random
from typing import Optional

import numpy as np

import ase
from ase.utils.ff import Morse, Angle, Dihedral, VdW
from ase.calculators.ff import ForceField
from ase.optimize.precon import FF, PreconLBFGS
from ase.optimize.precon.neighbors import get_neighbours

from rdkit import Chem
from rdkit.Chem.AllChem import (
    EmbedMultipleConfs, GetBestRMS,
    MMFFGetMoleculeForceField, MMFFGetMoleculeProperties,
    UFFGetMoleculeForceField,
)

from . import settings


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


def calc_ff(atoms: ase.Atoms) -> tuple[ForceField, FF]:
    neighbor_list = [[] for _ in range(len(atoms))]
    vdw_list = np.ones((len(atoms), len(atoms)), dtype=bool)
    morses = []; angles = []; dihedrals = []; vdws = []

    i_list, j_list, d_list, fixed_atoms = get_neighbours(atoms=atoms, r_cut=1.5)
    for i, j in zip(i_list, j_list):
        neighbor_list[i].append(j)
    for i in range(len(neighbor_list)):
        neighbor_list[i].sort()

    for i in range(len(atoms)):
        for jj in range(len(neighbor_list[i])):
            j = neighbor_list[i][jj]
            if j > i:
                morses.append(Morse(atomi=i, atomj=j, D=6.1322, alpha=1.8502, r0=1.4322))
            vdw_list[i, j] = vdw_list[j, i] = False
            for kk in range(jj+1, len(neighbor_list[i])):
                k = neighbor_list[i][kk]
                angles.append(Angle(atomi=j, atomj=i, atomk=k, k=10.0, a0=np.deg2rad(120.0), cos=True))
                vdw_list[j, k] = vdw_list[k, j] = False
                for ll in range(kk+1, len(neighbor_list[i])):
                    l = neighbor_list[i][ll]
                    dihedrals.append(Dihedral(atomi=j, atomj=i, atomk=k, atoml=l, k=0.346))

    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            if vdw_list[i, j]:
                vdws.append(VdW(atomi=i, atomj=j, epsilonij=0.0115, rminij=3.4681))

    return (
        ForceField(morses=morses, angles=angles, dihedrals=dihedrals, vdws=vdws),
        FF(morses=morses, angles=angles, dihedrals=dihedrals),
    )


def pre_optimize(atoms: ase.Atoms, fmax: float = 1e-4) -> None:
    logfile = Path(settings.SCRATCH_PATH).joinpath("temp.pre")

    calc, precon = calc_ff(atoms)

    #_calc = atoms.calc
    atoms.calc = calc
    opt = PreconLBFGS(atoms, precon=precon, use_armijo=True, logfile=logfile)
    opt.run(fmax=fmax)
    #atoms.calc = _calc
    return atoms

