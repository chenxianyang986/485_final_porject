from ase.build import bulk
from ase import Atoms
from ase.atom import Atom
from ase.calculators.emt import EMT
from ase.phonons import Phonons
from ase.dft.kpoints import monkhorst_pack
import matplotlib.pyplot as plt
import numpy as np
from ase.calculators.eam import EAM
from ase.md import verlet
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.calculators.lj import LennardJones
import meam_calculator as meamc
import eam_calculator as eamc
import sys
from copy import deepcopy
import math

h_bar = 6.582 * 10**-16
def G_bcc(T,V, calc):

    temp = str(T)
    vol = str(V)
    tempname = (vol+"_"+temp + "_phonon" +  "_bcc")
    
    a_ = (2*V)**(1/3)
    atoms = bulk('Zr', 'bcc',  a=a_)
    if calc == 'meam':
        atoms.calc = meamc.setup_meam_calculator(atoms)
    else:
        atoms.calc = eamc.set_up_eam_calculator()
        print("get here!")
    N = 7
    ph = Phonons(atoms, atoms.calc, supercell=(N, N, N), delta=0.05,name = tempname)
    ph.run()
    ph.read(acoustic=True)
    ph.clean()

    #Trial1:
    kv=(7, 7, 7)
    kvc = monkhorst_pack(kv)
    omega_kl = ph.band_structure(kvc,verbose = False)
    return omega_kl

def G_hcp(T,V, calc):
  
  a_ = (V/np.sqrt(2))**(1/3)
  
  temp = str(T)
  vol = str(V)
  tempname = (vol+"_"+temp + "_phonon" +  "_hcp")
  
  atoms = bulk('Zr', 'hcp',  a=a_, c = 1.633*a_)
  if calc == 'meam':
        atoms.calc = meamc.setup_meam_calculator(atoms)
  else:
        atoms.calc = eamc.set_up_eam_calculator()
  N = 7
  ph = Phonons(atoms, atoms.calc, supercell=(N, N, N), delta=0.05,name = tempname)
  ph.run()
  ph.read(acoustic=True)
  ph.clean()

  #Trial1:
  kv=(7, 7, 7)
  kvc = monkhorst_pack(kv)
  omega_kl = ph.band_structure(kvc,verbose = False)
  
  return omega_kl

def differential_hcp(V,T, calc):
    
    V1 = V+0.005
    V2 = V-0.005
    omega1 = np.array(G_hcp(T,V1, calc))
    omega2 = np.array(G_hcp(T,V2, calc))
    P = (omega1 - omega2)/(0.01 * h_bar)
    
    return P

def differential_bcc(V,T, calc):
    
    V1 = V+0.005
    V2 = V-0.005
    omega1 = np.array(G_bcc(T,V1, calc))
    omega2 = np.array(G_bcc(T,V2, calc))
    P = (omega1 - omega2)/(0.01 * h_bar)
    
    return P

def cvt(omega_kl,T):
    k = 8.617333262e-5
    w = np.abs(deepcopy(omega_kl))
    d = k*(w/(k*T))**2*(np.exp(w/(k*T))/(np.exp(w/(k*T)) - 1)**2) 
    return d

def Gruneisen_hcp(V,T, calc):
    omega_kl = G_hcp(T, V, calc)
    omega = np.array(omega_kl)
    diff = differential_hcp(V,T, calc)
    g = diff*(V/(omega/h_bar))
    c = cvt(omega,T)
    gamma = np.sum(np.multiply(g,c))/np.sum(c)
    return gamma,g

def Gruneisen_bcc(V,T, calc):
    
    omega_kl = G_bcc(T, V, calc)
    omega = np.array(omega_kl)
    diff = differential_bcc(V,T, calc)
    g = diff*V/(omega/h_bar)
    c = cvt(omega,T)
    gamma = np.sum(np.multiply(g,c))/np.sum(c)
    return gamma,g

def main(argv):
    T = [200, 400, 600, 800]
    lattice_constant_a = np.linspace(3.22, 3.26, len(T))
    lattice_constant_c = np.linspace(5.21, 5.28, len(T))
    gamma_bccs = []
    gamma_hcps = []
    for i in range(len(T)):
        atoms_bcc: Atoms = bulk('Zr', 'bcc', a=lattice_constant_a[i])
        atoms_hcp: Atoms = bulk('Zr', 'hcp', a=lattice_constant_a[i], c=lattice_constant_c[i])
        calc1 = "meam"
        if argv[0] == "--eam":
            calc1 = "eam"
        gamma_hcp, _ = Gruneisen_hcp(atoms_hcp.get_volume(), T[i], calc1)
        gamma_bcc, _ = Gruneisen_bcc(atoms_bcc.get_volume(), T[i], calc1)
        gamma_bccs.append(abs(gamma_bcc))
        gamma_hcps.append(abs(gamma_hcp))
    #print(gamma_bccs, gamma_hcps)
    
    plt.plot(T, gamma_hcps, label = 'hcp')
    plt.plot(T, gamma_bccs, label = 'bcc')
    plt.xlabel("Temperature K")
    plt.ylabel("Grunesein gamma value")
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    main(sys.argv[1:])