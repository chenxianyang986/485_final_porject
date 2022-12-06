#!/usr/bin/env python
# coding: utf-8

# In[51]:


#sinusoid
from ase.build import bulk
from ase import Atoms
from ase.atom import Atom
from ase.calculators.emt import EMT
from ase.phonons import Phonons
from ase.dft.kpoints import *
import matplotlib.pyplot as plt
import numpy as np
from ase.calculators.eam import EAM
from ase.md import verlet
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.calculators.lj import LennardJones


# In[49]:


def MC_2part(d,calc,steps = 600,stepsize = 1):
    a = Atoms([Atom('Zr', (0, 0, -d/2)), Atom('Zr', (0, 0, d/2))],calculator = calc)

    dyn = VelocityVerlet(a,timestep=stepsize*ase.units.fs)

    forces = []
    Etotal = []

    def forceplot(a):
        F = a.get_forces() 
        forces.append(F[0][2])
        etot = a.get_potential_energy() + a.get_kinetic_energy()
        Etotal.append(etot)


    for i in range(steps):
        dyn.run(1)
        forceplot(a)
    
    return forces, Etotal


def periodcheck(f):
    ps = []
    m1 = np.min(np.array(f))
    
    for i in  range(len(f)):
        if(f[i]-m1 < 0.5):
            ps.append(i)
            
    return ps

def periodiso(f,ps):
    for i in range(len(ps-1)):
        if((ps[i+1] - ps[i]) > 20):
            return f[ps[i],ps[i+1]]
        
    return f


# In[50]:


def HAM_check(d,calc,steps = 600,stepsize = 1):
    
    forces,energy = MC_2part(d,calc,steps = 600,stepsize = 1)
    ps = periodcheck(forces)
    f = periodiso(forces,ps)
    
    dft = np.abs(np.fft.rfft(f))
    ndft = np.abs(dft)/np.trapz(np.abs(dft))
    
    Ham = np.std(ndft)**2
    
    return forces,f,dft,ndft,Ham


# In[ ]:





# In[ ]:





# In[ ]:




