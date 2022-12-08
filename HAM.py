from ase.build import bulk
from ase import Atoms
from ase.atom import Atom
from ase.calculators.emt import EMT
from ase.phonons import Phonons
from ase.dft.kpoints import monkhorst_pack
import matplotlib.pyplot as plt
from ase.visualize import view
import numpy as np
import eam_calculator as eamc
from ase.md import verlet
from ase.md.verlet import VelocityVerlet
from ase import units
import meam_calculator as meamc
import sys

def periodcheck(f):
    ps = []
    m1 = np.min(np.array(f))
    
    for i in  range(len(f)):
        if(f[i]-m1 < 0.3):
            ps.append(i)
            
    return ps

def periodiso(f,ps):
    for i in range(len(ps)-1):
        if((ps[i+1] - ps[i]) > 20):
            return f[ps[i]:ps[i+1]]
        
    return f

def advance(atoms,calc,idx,vel,dt):
    #initialize
    atoms.calc = calc
    mass = np.array(atoms.get_masses())
    m = mass[3]
    force = np.array(atoms.get_forces())
    pos = np.array(atoms.get_positions())
    
    #move1
    accel = force[idx] / m
    vel_half = vel + 0.5*dt*accel
    pos_new = pos 
    pos_new[idx]+= dt*vel_half
    
    atoms.set_positions(pos_new)
    
    #move2
    force2 = np.array(atoms.get_forces())
    accel2 = force2[idx] / m
    vel_new = vel_half + 0.5*dt*accel2
    
    return force2[idx],vel_new,pos_new[idx]

def changearray_hcp(mode,dist,atom):
    c = np.zeros(shape = (54,3))
    c[atom] = (mode)*dist
    return c

def GridMC_hcp_1(calc,steps,mode,dist,atom):
    dt = 0.1
    
    atoms = bulk('Zr', 'hcp', a=3.24, c = 5.22) * (3, 3, 3)
    j = np.array(atoms.get_positions())
    cinit = changearray_hcp(mode,dist,atom)
    pos1 = j + cinit
    atoms.set_positions(pos1)
    
    vel = [0,0,0]
    
    forcearray = []
    posarray = []
    
    for i in range(steps):
        f,vel,pos = advance(atoms,calc,atom,vel,dt)
        forcearray.append(f)
        posarray.append(pos)
        if(i%20 ==0):
            print(i)
        
    return forcearray,posarray


def GridMC_bcc_1(calc,steps,mode,dist,atom):
    dt = 0.1
    
    atoms = bulk('Zr', 'bcc', a=3.6) * (3, 3, 3)
    j = np.array(atoms.get_positions())
    cinit = changearray_hcp(mode,dist,atom)
    pos1 = j + cinit
    atoms.set_positions(pos1)
    
    vel = [0,0,0]
    
    forcearray = []
    posarray = []
    
    for i in range(steps):
        f,vel,pos = advance(atoms,calc,atom,vel,dt)
        forcearray.append(f)
        posarray.append(pos)
        if(i%20 ==0):
            print(i)
        
    return forcearray,posarray

def HAM_check_multi_hcp(calc, displacement=1,mode=[1,1,0],atom=14,steps=600):
    
    forc,posi = GridMC_hcp_1(calc,steps,mode,displacement,atom)
    farray = np.array(forc)
    forces = np.sqrt(farray[:,0]**2 + farray[:,1]**2+farray[:,2]**2)
    
    ps = periodcheck(forces)
    f = periodiso(forces,ps)
    
    dft = np.abs(np.fft.rfft(f))
    ndft = np.abs(dft)/np.trapz(np.abs(dft))
    
    Ham = np.std(ndft)**2
    
    return forces,f,dft,ndft,Ham
    
    
    
def HAM_check_multi_bcc(calc, displacement=1,mode=[1,1,0],atom=14,steps=600):
    
    forc,posi = GridMC_bcc_1(calc,steps,mode,displacement,atom)
    farray = np.array(forc)
    forces = np.sqrt(farray[:,0]**2 + farray[:,1]**2+farray[:,2]**2)
    
    ps = periodcheck(forces)
    f = periodiso(forces,ps)
    
    dft = np.abs(np.fft.rfft(f))
    ndft = np.abs(dft)/np.trapz(np.abs(dft))
    
    Ham = np.std(ndft)**2
    
    return forces,f,dft,ndft,Ham

def main(argv):
  if argv[0] == "--meam_single":
    Calc = meamc.setup_meam_calculator(None)
    f,p = GridMC_hcp_1(Calc,600,[1,1,0],1,14)
    plt.plot(f)
    plt.xlabel("steps")
    plt.ylabel("Force")
    plt.show()
  elif argv[0] == "--all_potentials":
    forces1,f1,dft1,ndft1,Ham1 = HAM_check_multi_hcp(displacement=1,mode=[1,1,0],atom=14,steps=600,calc=eamc.set_up_leonard_jones_calculator())
    forces2,f2,dft2,ndft2,Ham2 = HAM_check_multi_hcp(displacement=1,mode=[1,1,0],atom=14,steps=600,calc=eamc.set_up_eam_calculator())
    forces4,f4,dft4,ndft4,Ham4 = HAM_check_multi_hcp(displacement=1,mode=[1,1,0],atom=14,steps=600,calc=meamc.setup_meam_calculator(None))
    print(Ham1, Ham2, Ham4)
    '''
    plt.plot(forces1, label="LJ")
    plt.plot(forces2, label="EAM")
    plt.plot(forces4, label="MEAM")
    plt.legend()
    plt.xlabel("steps")
    plt.ylabel("Force")
    plt.show()
    '''
    plt.plot(ndft1, label="LJ")
    plt.plot(ndft2, label="EAM")
    plt.plot(ndft4, label="MEAM")
    plt.legend()
    plt.xlabel("steps")
    plt.ylabel("normalized_Fourier_transformed_force")
    plt.show()
if __name__ == "__main__":
  main(sys.argv[1:])