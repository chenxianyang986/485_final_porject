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
number_of_sweeps = 300
vmax = 2
beta = 1
pressure = 1
mu = 0
tau = 1
calc = eamc.set_up_eam_calculator()
atoms = ba.set_bcc_atoms_in_volume(lengths, int(n_cells ** (1/3)))
n_atoms = len(atoms)
atoms.calc = calc
potential = []
volume_acceptance_time = 0
total_volume_attempt = 0
position_acceptance_ratio = []
all_atom_positions = []
for i in range(number_of_sweeps):
    print(i)
    moves = np.random.normal(mu, tau, n_atoms * 3)
    moves = moves.reshape((n_atoms, 3))
    acceptance_check = np.random.uniform(size = n_atoms)
    volume_accept, move_accept, current_potential, state, positions = mc.my_mc_sweep(atoms, lengths, beta, moves, acceptance_check, vmax, pressure)
    all_atom_positions.append(positions)
    if state == True:
        total_volume_attempt += 1
        if volume_accept:
            volume_acceptance_time += 1
    else:
        position_acceptance_ratio.append(move_accept)
    potential.append(current_potential)
if total_volume_attempt != 0:
    print("volume update acc is", volume_acceptance_time / total_volume_attempt)
else:
    print("no volume change occur in this time period")
pcf.get_pair_correlation(1/beta, all_atom_positions, n_atoms, lengths)



