from ase.calculators.lj import LennardJones
import build_atom as ba
from ase.visualize import view
import eam_calculator as eamc
import matplotlib.pyplot as plt

calc = eamc.set_up_eam_calculator()
lengths = 10

atoms = ba.set_bcc_atoms_in_volume(lengths, 3)
#atoms.set_pbc(True)
atoms.calc = calc
atoms.get_potential_energy()


