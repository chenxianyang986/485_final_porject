import meam_calculator as meamc
import numpy as np
from ase import Atoms, Atom
import matplotlib.pyplot as plt
import sys
from ase.build import bulk
from ase.visualize import view

def test_meam_potential():
  calc = meamc.setup_meam_calculator()
  lengths = np.linspace(1, 5, 301)
  potentials = []
  atomic_distances = []

  
  for i in range(len(lengths)):
    atoms: Atoms = bulk('Zr', crystalstructure='hcp')
    atoms.set_cell((atoms.get_cell().lengths() * i).tolist() + atoms.get_cell().angles().tolist(), scale_atoms=True)
    print(atoms.get_cell())
    atoms.calc = calc
    potentials.append(atoms.get_potential_energy())
    atomic_distances.append(lengths[i])
  plt.plot(atomic_distances, potentials, color='green', linestyle='dashed')
  plt.xlabel("r Angstrom")
  plt.ylabel("Energy ev")
  plt.show()
  
  #Ni = bulk('Ni', cubic=True)
  #Fe = Atom('Fe', position=Ni.cell.diagonal()/2)
  #NiFe = Ni + Fe
  #NiFe.calc = calc
  #print("Energy ", NiFe.get_potential_energy())
  #Zr_1 = Atom('Fe', position=Zr.cell.diagonal()/2)
  #Zr2 = Zr + Zr_1
  #Zr2.calc = calc
  #print("Energy ", Zr2.get_potential_energy())

def main(argv):
  if argv[0] == "--test_meam":
    test_meam_potential()

if __name__ == "__main__":
  main(sys.argv[1:])