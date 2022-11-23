import meam_calculator as meamc
import numpy as np
from ase import Atoms
import matplotlib.pyplot as plt
import sys

def test_meam_potential():
  calc = meamc.setup_meam_calculator()
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

def main(argv):
  if argv[0] == "--test_meam":
    test_meam_potential()

if __name__ == "__main__":
  main(sys.argv[1:])