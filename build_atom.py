from ase.build import bulk
from ase.visualize import view
from ase import Atoms
import numpy as np
from ase.calculators.lj import LennardJones
import math

def set_bcc_atoms_in_volume(length, cell_number_one_side):

    a1: Atoms = bulk('Zr', 'bcc', a = length / cell_number_one_side, cubic=True)
    basic_atom_positions = np.array(a1.get_positions())
    one_side = np.linspace(-length/2, length/2, cell_number_one_side + 1)[:-1]
    x, y, z = np.meshgrid(one_side, one_side, one_side)
    atom_pos_in_box = []
    for i in range(len(x.flatten())):
        for j in range(len(basic_atom_positions)):
            atom_pos_in_box.append((np.array([x.flatten()[i], y.flatten()[i], z.flatten()[i]]) + np.array(basic_atom_positions[j])))
    cell_parameter = [length, length, length]
    atoms = Atoms(len(atom_pos_in_box) * 'Zr', atom_pos_in_box)
    return atoms

def set_hcp_atoms_in_volume(base_length, height_length, _a, _c):
    a1: Atoms = bulk('Zr', 'hcp', a = base_length / _a, c = height_length / _c)
    basic_atom_positions = np.array(a1.get_positions())
    base_side = np.linspace(-base_length / 2, base_length/2, int(_a) + 1)[:-1]
    height_side = np.linspace(-height_length / 2, height_length/2, int(_c) + 1)[:-1]
    x, y, z = np.meshgrid(height_side, base_side, base_side)
    atom_pos_in_box = []
    for i in range(len(x.flatten())):
        if ((i // int(_a)) % int(_a)) % 2 == 1:
            for k in range(len(basic_atom_positions)):
                atom_pos_in_box.append((np.array([x.flatten()[i], y.flatten()[i] + height_length / int(_c) / 2, z.flatten()[i] + height_length / int(_c) / 2]) + np.array(basic_atom_positions[k])))
        else:
            for k in range(len(basic_atom_positions)):
                atom_pos_in_box.append((np.array([x.flatten()[i], y.flatten()[i], z.flatten()[i]]) + np.array(basic_atom_positions[k])))
    atoms = Atoms('Zr' * len(atom_pos_in_box), positions=atom_pos_in_box)
    return atoms