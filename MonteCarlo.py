import math
from ase import Atoms
import eam_calculator as eamc
import numpy as np
import math 
from copy import deepcopy

def check_distance(cur_pos, all_pos, index):
  for i in range(len(all_pos)):
    distance = np.sum((cur_pos - all_pos[i]) ** 2)
    if distance < 0.01 and index != i:
      return False
  return True

def my_loga_symmetric(vold, vnew, beta):
  """ Calculate the logarithm of the acceptance probability
   given the energies before and after a symmetric move.

  Args:
    vold (float): potential energy before move
    vnew (float): potential energy after move
    beta (float): inverse temperature
  Return:
    float: ln(A(x->x')), logorithm of the acceptance probability
  """
  return - beta * (vnew - vold)

def volume_change_acceptance_rule(potential_old, potential_new, volume_old, volume_new, beta, pressure, particle_number):
  return particle_number * np.log(volume_new / volume_old) - beta * (potential_new - potential_old) - beta * pressure * (volume_new - volume_old)

def my_mc_sweep(atoms: Atoms, lbox, beta, eta, acc_check, vmax, pressure):
  """ Perform one Monte Carlo sweep

  Args:
    atoms: ASE atoms with a calc 
    lbox (float): cubic box side length
    beta (float): inverse temperature
    eta (np.array): shape (natom, ndim), array of Gaussian random
     numbers, one for each single-particle move
    acc_check (np.array): shape (natom), array of uniform random
     numbers, one for each single-particle move
    vmax(float): maximum to update volume
    pressure(float): a fixed value 
  Return:
    (int, float): (naccept, de), the number of accepted
    single-particle moves, and the change in energy
  """
  volume = lbox ** 3
  alpha = np.random.random()
  pos = atoms.get_positions()
  natom, ndim = pos.shape
  n_volume_accept = False
  n_move_accept = 0
  # complete this function
  oldpos = deepcopy(pos)
  state = False # False for a position move; otherwise, true
  if alpha < 1/(natom + 1):
    state = True
    oldpotential = atoms.get_potential_energy()
    new_volume = math.e **((np.random.random() - 0.5) * vmax + np.log(volume))
    new_lbox = new_volume ** (1/3)
    new_pos = deepcopy(oldpos) * new_lbox / lbox
    atoms.set_positions(eamc.minimum_image(new_pos, new_lbox))
    newpotential = atoms.get_potential_energy()
    acc_probability = math.e ** volume_change_acceptance_rule(oldpotential, newpotential, volume, new_volume, beta, pressure, natom)
    if np.random.random() < acc_probability:
      n_volume_accept = True
    else:
      atoms.set_positions(eamc.minimum_image(oldpos, lbox))
  else:
    for i in range(natom):
      oldpotential = atoms.get_potential_energy()
      newpos = deepcopy(oldpos)
      for j in range(ndim):
        newpos[i][j] += eta[i][j]
      atoms.set_positions(eamc.minimum_image(newpos, lbox))
      newpotential = atoms.get_potential_energy()
      acc_probability = math.e ** my_loga_symmetric(oldpotential, newpotential, beta)
      if acc_check[i] <= acc_probability and newpotential > -2000: #TODO: may change later for this check
        n_move_accept += 1
        oldpos = deepcopy(eamc.minimum_image(newpos, lbox))
      else:
        atoms.set_positions(eamc.minimum_image(oldpos, lbox))
  return n_volume_accept, n_move_accept / natom, atoms.get_potential_energy(), state, atoms.get_positions()
