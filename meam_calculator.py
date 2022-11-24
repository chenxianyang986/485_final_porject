from ase.calculators.lammpsrun import LAMMPS
import os
from ase.calculators.lammpslib import LAMMPSlib

import lammps
def setup_meam_calculator() -> LAMMPSlib:
  '''
  parameters = {'pair_style': 'meam', 'pair_coeff': ['* * library.meam.fs Zr']}
  files = {"Zr.meam.fs"}
  return LAMMPS(parameters = parameters, files=files)
  '''
  cmd = ["pair_style meam",
        "pair_coeff * * library.meam Zr NULL Zr"]
  return LAMMPSlib(lmpcmds = cmd)