from ase.calculators.lj import LennardJones
import build_atom as ba
from ase.visualize import view
import eam_calculator as eamc
import matplotlib.pyplot as plt
import sys
import numpy as np
import MonteCarlo as mc

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
  for i in range(number_of_sweeps):
    #print(i)
    moves = np.random.normal(mu, tau, n_atoms * 3)
    moves = moves.reshape((n_atoms, 3))
    acceptance_check = np.random.uniform(size = n_atoms)
    volume_accept, move_accept, current_potential, state = mc.my_mc_sweep(atoms, lengths, beta, moves, acceptance_check, vmax, pressure)
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
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
  axes[0].plot(position_acceptance_ratio)
  axes[1].plot(potential)
  fig.tight_layout()
  plt.show()

def main(argv):
  if argv[0] == "--setup":
    test_atom_set_up()
  if argv[0] == "--pbc":
    test_intial_set_up()
  if argv[0] == "--mc":
    test_monte_carlo()

if __name__ == "__main__":
    main(sys.argv[1:])