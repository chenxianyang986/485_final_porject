from ase.calculators.eam import EAM
import numpy as np
def set_up_eam_calculator() -> EAM:
    calculator = EAM(potential = "./Zr_2.eam.fs")
    return calculator

def minimum_image(r, L):
  '''
  r position of atom
  L length of the box
  return position of atom in box
  '''
  return r - L*np.round(r / L)