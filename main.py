from ase.calculators.lj import LennardJones
import build_atom as ba
from ase.visualize import view
import eam_calculator as eamc
import matplotlib.pyplot as plt
#atoms = ba.set_bcc_atoms_in_volume(7.2, 3.6)

calc = eamc.set_up_eam_calculator()
#LennardJones(sigma = 1, epsilon = 1, rc = 7.2 /1.5)
#atoms.calc = calc
#print(atoms.get_potential_energy())
#print(atoms.get_volume() / atoms.get_global_number_of_atoms())
#view(atoms)
lengths = [4.5, 5, 5.6, 6, 7.2, 9.6, 12, 14.4, 16, 17, 18, 20, 24]
potentials = []
atomic_volumes = []
for i in range(len(lengths)):
    atoms = ba.set_bcc_atoms_in_volume(lengths[i], lengths[i] / 2)
    atoms.calc = calc
    potentials.append(atoms.get_potential_energy())
    atomic_volumes.append(atoms.get_volume() / atoms.get_global_number_of_atoms())
    #view(atoms)
plt.plot(atomic_volumes, potentials, color='green', marker='o', linestyle='dashed')
plt.show()
'''
a1 = ba.set_hcp_atoms_in_volume(2.4, 7.2, 1.2, 3.6)
view(a1)
'''