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
import structure_factor as sf

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
  lengths = np.linspace(2, 8, 501)
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
  lengths = 3.24 * 3
  n_cells = 27
  number_of_sweeps = 2500
  vmax = 1.8
  beta = 1 / (8.617 * 10 ** -5 * 750)
  pressure = 0.01
  mu = 0
  tau = 0.65
  calc = eamc.set_up_eam_calculator()
  atoms = ba.set_bcc_atoms_in_volume(lengths, int(n_cells ** (1/3)))
  view(atoms)
  n_atoms = len(atoms)
  atoms.calc = calc
  potential = []
  volume_acceptance_time = 0
  total_volume_attempt = 0
  position_acceptance_ratio = []
  reduced_kvecs, sk_list = plot_sk(4, lengths, atoms.get_positions())
  plt.plot(reduced_kvecs[1:], sk_list[1:])
  plt.show()
  positions_total = []
  for i in range(number_of_sweeps):
    print(i)
    moves = np.random.normal(mu, tau, n_atoms * 3)
    moves = moves.reshape((n_atoms, 3))
    acceptance_check = np.random.uniform(size = n_atoms)
    volume_accept, move_accept, current_potential, state, positions = mc.my_mc_sweep(atoms, lengths, beta, moves, acceptance_check, vmax, pressure)
    #print(current_potential, state)
    if i >= number_of_sweeps * 2/3:
      positions_total.append(np.array(positions))
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
  total_sk = np.array([])
  reduced_kvecs = np.array([])
  for i in positions_total:
    reduced_kvecs, reduced_sk = plot_sk(4, lengths, i)
    if len(total_sk) == 0:
      total_sk = reduced_sk
    else:
      total_sk += reduced_sk
  plt.plot(reduced_kvecs[1:], (total_sk / len(positions_total))[1:])
  plt.xlabel("kvecs")
  plt.ylabel("sk")
  plt.show()

def test_monte_carlo_in_hcp():
  lengths = 3.24 * 3
  n_cells = 27
  number_of_sweeps = 30000 #must be larger than 500
  vmax = 0.5
  beta = 1/(8.617 * 10 ** -5 * 1100)
  pressure = 0.01
  mu = 0
  tau = 0.55
  calc = eamc.set_up_eam_calculator()
  atoms = ba.set_hcp_atoms_in_volume(lengths, lengths, np.round(n_cells ** (1/3)), np.round(n_cells ** (1/3)))
  view(atoms)
  n_atoms = len(atoms)
  '''
  gr = pcf.get_pair_correlation(1/beta, [atoms.get_positions()], n_atoms, lengths)
  plt.plot(gr)
  plt.show()
  '''
  reduced_sk, sk_list = plot_sk(4, lengths, atoms.get_positions())
  plt.plot(reduced_sk[1:], sk_list[1:])
  plt.xlabel("kvecs")
  plt.ylabel("sk")
  plt.show()
  atoms.calc = calc
  potential = []
  volume_acceptance_time = 0
  total_volume_attempt = 0
  position_acceptance_ratio = []
  positions_total = []

  for i in range(number_of_sweeps):
    print(i)
    moves = np.random.normal(mu, tau, n_atoms * 3)
    moves = moves.reshape((n_atoms, 3))
    acceptance_check = np.random.uniform(size = n_atoms)
    volume_accept, move_accept, current_potential, state, positions = mc.my_mc_sweep(atoms, lengths, beta, moves, acceptance_check, vmax, pressure)
    #print(current_potential, state)
    
    if i >= 200:
      positions_total.append(np.array(positions))
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
  '''
  gr_simulation = pcf.get_pair_correlation(1/beta, positions_total, n_atoms, lengths)
  plt.plot(gr_simulation)
  plt.show()
  '''
  total_sk = np.array([])
  reduced_kvecs = np.array([])
  for i in positions_total:
    reduced_kvecs, reduced_sk = plot_sk(4, lengths, i)
    if len(total_sk) == 0:
      total_sk = reduced_sk
    else:
      total_sk += reduced_sk
  plt.plot(reduced_kvecs[1:], (total_sk / len(positions_total))[1:])
  plt.xlabel("kvecs")
  plt.ylabel("sk")
  plt.show()

def plot_sk(maxn, box_length, avg_pos):
  kvecs = sf.my_legal_kvecs(maxn, box_length)
  sk = sf.my_calc_sk(kvecs, avg_pos)
  reduced_kvecs, reduced_sk = sf.structure_factor_plot_helper(kvecs, sk)
  return reduced_kvecs, reduced_sk

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
    plt.plot(gr_simulation, label= f"1/kT is {beta}")
    plt.xlabel("r*")
    plt.ylabel("gr")
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