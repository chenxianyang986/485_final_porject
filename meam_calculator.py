from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.lammpslib import LAMMPSlib
from lammps import lammps
lmp = lammps()

def setup_meam_calculator(atoms):
  cmd = ["pair_style meam",
        "pair_coeff * * library.meam Zr Zr.meam Zr"]
  return LAMMPSlib(lmpcmds = cmd, atoms=atoms, keep_alive=True)