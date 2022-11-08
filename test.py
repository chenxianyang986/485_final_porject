from ase.calculators.lj import LennardJones
import build_atom as ba
from ase.visualize import view
import eam_calculator as eamc
import matplotlib.pyplot as plt
import sys
import numpy as np
import MonteCarlo as mc
import pair_correlation_function as pcf
from ase import Atoms
from copy import deepcopy

np.set_printoptions(threshold=sys.maxsize)

def test_atom_set_up():
  calc = eamc.set_up_eam_calculator()
  lengths = [7.2, 8.4, 9.6, 12, 14.4, 16, 17, 18, 20, 24]
  potentials = []
  atomic_volumes = []
  for i in range(len(lengths)):
      atoms = ba.set_bcc_atoms_in_volume(lengths[i], 3)
      atoms.calc = calc
      print(atoms.positions, lengths[i])
      potentials.append(atoms.get_potential_energy())
      atomic_volumes.append(atoms.get_volume() / atoms.get_global_number_of_atoms())
      #view(atoms)
  plt.plot(atomic_volumes, potentials, color='green', marker='o', linestyle='dashed')
  plt.show()

def test_atom_set_up_2():
  calc = eamc.set_up_eam_calculator()
  lengths = np.linspace(1.8, 5, 301)
  potentials = []
  atomic_distances = []
  for i in range(len(lengths)):
      atoms = Atoms(2 * "Zr", positions=[[0, 0, 0], [0, 0, lengths[i]]])
      atoms.calc = calc
      potentials.append(atoms.get_potential_energy())
      atomic_distances.append(lengths[i])
      print(lengths[i], potentials[i])
      #view(atoms)
  plt.plot(atomic_distances, potentials, color='green', linestyle='dashed')
  plt.xlabel("r Angstrom")
  plt.ylabel("Energy ev")
  plt.show()

def test_intial_set_up():
  lengths = 20
  calc = eamc.set_up_eam_calculator()
  atoms = ba.set_bcc_atoms_in_volume(lengths, 4)
  atoms.calc = calc
  atoms_test = ba.set_bcc_atoms_in_volume(lengths, 4)
  atoms_test.calc = calc
  positions = atoms.get_positions()
  for i in range(len(positions)):
      positions[i] = eamc.minimum_image(np.array(positions[i]), lengths)
  test_positions = atoms_test.get_positions()
  print(positions)
  assert len(positions) == len(test_positions)
  for i in range(len(positions)):
      assert all(positions[i] == test_positions[i])
  print("initial set up for pbc is okay!")
  print(atoms.get_potential_energy())

def test_monte_carlo():
  lengths = 20
  n_cells = 27
  number_of_sweeps = 500
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
  for i in range(number_of_sweeps):
    print(i)
    moves = np.random.normal(mu, tau, n_atoms * 3)
    moves = moves.reshape((n_atoms, 3))
    acceptance_check = np.random.uniform(size = n_atoms)
    volume_accept, move_accept, current_potential, state, positions = mc.my_mc_sweep(atoms, lengths, beta, moves, acceptance_check, vmax, pressure)
    #print(current_potential, state)
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
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
  axes[0].plot(position_acceptance_ratio)
  axes[0].set_xlabel("number of sweeps")
  axes[0].set_ylabel("acceptance ratio")
  axes[1].plot(potential)
  axes[1].set_xlabel("number of sweeps")
  axes[1].set_ylabel("potential")
  fig.tight_layout()
  plt.show()
  view(atoms)

def test_monte_carlo_in_hcp():
  lengths = 3.24 * 3
  n_cells = 27
  number_of_sweeps = 200
  vmax = 2
  beta = 1/(8.617 * 10 ** -5 *1000)
  pressure = 1 * 10** 5
  mu = 0
  tau = 0.3
  calc = eamc.set_up_eam_calculator()
  atoms = ba.set_hcp_atoms_in_volume(lengths, lengths, int(n_cells ** (1/3)), int(n_cells ** (1/3)))
  view(atoms)
  n_atoms = len(atoms)
  atoms.calc = calc
  potential = []
  volume_acceptance_time = 0
  total_volume_attempt = 0
  position_acceptance_ratio = []
  #positions_total = np.array([])
  for i in range(number_of_sweeps):
    print(i)
    moves = np.random.normal(mu, tau, n_atoms * 3)
    moves = moves.reshape((n_atoms, 3))
    acceptance_check = np.random.uniform(size = n_atoms)
    volume_accept, move_accept, current_potential, state, positions = mc.my_mc_sweep(atoms, lengths, beta, moves, acceptance_check, vmax, pressure)
    #print(current_potential, state)
    '''
    if len(positions_total) == 0 and i >= number_of_sweeps - 50:
      positions_total = deepcopy(positions)
    elif  i >= number_of_sweeps - 50 and len(positions_total) != 0:
      positions_total += positions
    '''
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
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
  axes[0].plot(position_acceptance_ratio)
  axes[0].set_xlabel("number of sweeps")
  axes[0].set_ylabel("acceptance ratio")
  axes[1].plot(potential)
  axes[1].set_xlabel("number of sweeps")
  axes[1].set_ylabel("potential")
  fig.tight_layout()
  plt.show()
  #atoms.set_positions(positions_total / (50))
  view(atoms)

def test_monte_carlo_in_lenard_Jones():
  lengths = 5
  n_cells = 27
  number_of_sweeps = 100
  vmax = 2
  beta = 2
  pressure = 1
  mu = 0
  tau = 0.08
  #calc = eamc.set_up_eam_calculator()
  calc = eamc.set_up_leonard_jones_calculator()
  atoms = ba.set_bcc_atoms_in_volume(lengths, int(n_cells ** (1/3)))
  n_atoms = len(atoms)
  atoms.calc = calc
  potential = []
  volume_acceptance_time = 0
  total_volume_attempt = 0
  position_acceptance_ratio = []
  for i in range(number_of_sweeps):
    #print(i)
    moves = np.random.normal(mu, tau, n_atoms * 3)
    moves = moves.reshape((n_atoms, 3))
    acceptance_check = np.random.uniform(size = n_atoms)
    volume_accept, move_accept, current_potential, state, position= mc.my_mc_sweep(atoms, lengths, beta, moves, acceptance_check, vmax, pressure)
    print(current_potential, state)
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
  view(atoms)
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
  axes[0].plot(position_acceptance_ratio)
  axes[1].plot(potential)
  fig.tight_layout()
  plt.show()

def test_at_different_temperatures(betas):
  for beta in betas:
    gr_simulation = gr_simulation_helper_function(float(beta))
    plt.plot(gr_simulation, label="beta = {}".format(beta))
  plt.show()
  pass

def gr_simulation_helper_function(b):
  lengths = 20
  n_cells = 27
  number_of_sweeps = 300
  vmax = 2
  beta = b
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
  avg_pair_correlation = pcf.get_pair_correlation(1/beta, all_atom_positions, n_atoms, lengths)
  return avg_pair_correlation

def main(argv):
  if argv[0] == "--setup":
    test_atom_set_up()
  elif argv[0] == "--setup2":
    test_atom_set_up_2()
  elif argv[0] == "--pbc":
    test_intial_set_up()
  elif argv[0] == "--mc":
    test_monte_carlo()
  elif argv[0] == "--mchcp":
    test_monte_carlo_in_hcp()
  elif argv[0] == "--lj":
    test_monte_carlo_in_lenard_Jones()
  if argv[0]== "--plotgr":
    if len(argv) == 1:
      print("please provide your simulation temps!")
    else:
      test_at_different_temperatures(argv[1:])
if __name__ == "__main__":
    main(sys.argv[1:])