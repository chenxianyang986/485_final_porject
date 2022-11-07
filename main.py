from ase.calculators.lj import LennardJones
import build_atom as ba
from ase.visualize import view
import eam_calculator as eamc
import matplotlib.pyplot as plt
import numpy as np
import MonteCarlo as mc
import pair_correlation_function as pcf

lengths = 20
n_cells = 27

calc = eamc.set_up_eam_calculator()
atoms = ba.set_bcc_atoms_in_volume(lengths, int(n_cells ** (1/3)))
n_atoms = len(atoms)
atoms.calc = calc
view(atoms)