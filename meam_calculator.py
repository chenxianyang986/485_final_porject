from ase.calculators.lammpsrun import LAMMPS
import os
from ase.calculators.lammpslib import LAMMPSlib

def setup_meam_calculator() -> LAMMPSlib:
  '''
  os.environ["LAMMPSPATH"] = "C:\\Users\\xianyang chen\\mylammps"
  imp_path = os.environ['LAMMPSPATH']
  public_file_name = os.path.join(imp_path, "library.meam")  # type: ignore
  '''
  '''
  parameters = {'pair_style': 'meam', 'pair_coeff': ['* * library.meam.fs Zr']}
  #local_file_name = os.path.join(imp_path, "Zr.meam")  # type: ignore
  files = {"Zr.meam.fs"}
  return LAMMPS(parameters = parameters, files=files)
  '''
  cmd = {'pair_style': "lj/cut 2.5", 'pair_coeff':["* * 1 1"]}
  return LAMMPSlib(lmpcmd = cmd)