from ase.calculators.lj import LennardJones
import build_atom as ba
from ase.visualize import view

atoms = ba.set_bcc_atoms_in_volume(7.2, 3.6)

calc = LennardJones(sigma = 1, epsilon = 1, rc = 7.2 /1.5)
atoms.calc = calc
print(atoms.get_potential_energy())
view(atoms)

a1 = ba.set_hcp_atoms_in_volume(2.4, 7.2, 1.2, 3.6)
view(a1)