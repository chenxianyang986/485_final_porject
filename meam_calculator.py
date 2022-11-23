from ase.calculators.lammpsrun import LAMMPS
import os
def setup_meam_calculator() -> LAMMPS:
  '''
  os.environ["LAMMPSPATH"] = "C:\\Users\\xianyang chen\\mylammps"
  imp_path = os.environ['LAMMPSPATH']
  public_file_name = os.path.join(imp_path, "library.meam")  # type: ignore
  '''
  parameters = {'pair_style': 'meam', 'pair_coeff': ['* * library.meam Zr']}
  #local_file_name = os.path.join(imp_path, "Zr.meam")  # type: ignore
  files = {"Zr.meam"}
  return LAMMPS(parameters = parameters, files=files)