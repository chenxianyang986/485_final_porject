from ase.build import bulk
from ase.visualize import view
from ase import Atoms
import numpy as np
from ase.calculators.lj import LennardJones

def set_bcc_atoms_in_volume(length, lattice_constant):

    a1: Atoms = bulk('Zr', 'bcc', a = lattice_constant, cubic=True)
    basic_atom_positions = np.array(a1.get_positions())
    one_side = np.linspace(0, length, int(length / lattice_constant) + 1)
    x, y, z = np.meshgrid(one_side, one_side, one_side)
    atom_pos_in_box = []
    for i in range(len(x.flatten())):
        for j in range(len(basic_atom_positions)):
            atom_pos_in_box.append((np.array([x.flatten()[i], y.flatten()[i], z.flatten()[i]]) + np.array(basic_atom_positions[j])))
    atoms = Atoms(len(atom_pos_in_box) * 'Zr', atom_pos_in_box)
    return atoms

def set_hcp_atoms_in_volume(base_length, height_length, _a, _c):
    a1: Atoms = bulk('Zr', 'hcp', a = _a, c = _c)
    basic_atom_positions = np.array(a1.get_positions())
    base_side = np.linspace(0, base_length, int(base_length / _a) + 1)
    height_side = np.linspace(0, height_length, int(height_length / _c) + 1)
    x, y, z = np.meshgrid(base_side, base_side, height_side)
    atom_pos_in_box = []
    for i in range(len(x.flatten())):
            for k in range(len(basic_atom_positions)):
                atom_pos_in_box.append((np.array([x.flatten()[i], y.flatten()[i], z.flatten()[i]]) + np.array(basic_atom_positions[k])))
    atoms = Atoms(len(atom_pos_in_box) * 'Zr', atom_pos_in_box)
    return atoms