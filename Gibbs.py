from ase.build import bulk
from ase.calculators.emt import EMT
from ase.phonons import Phonons
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.dft.kpoints import monkhorst_pack
from ase.visualize import view
import math
import eam_calculator as eamc

h = 6.5821e-16
k_boltz = 8.617333262e-5


#METHOD 1: Discrete:
def gibbs(atoms,calc,N,P,V,T,kvecs = (10, 10, 10)):

    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    ph.read()
    ph.clean()

    Helm = 0
    #path = atoms.cell.bandpath('GXULGK', npoints=100)
    k = monkhorst_pack(kvecs)
    #print(path.kpts)

    bs = ph.band_structure(k)

    #print(bs)

    for i in range(len(bs)):
        for j in range(len(bs[i])):
            if(bs[i][j] > 0):
                Helm+= k*T*np.log(2*np.sinh(bs[i][j]/(2*k*T)))

    return(Helm/(N**3) + P*V + atoms.get_potential_energy())

def toy_gibbs(atoms,calc,N,P,V,T,kvecs = (6, 6, 6)):

    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    ph.read()
    ph.clean()

    Helm = 0
    #path = atoms.cell.bandpath('GXULGK', npoints=100)
    k = monkhorst_pack(kvecs)
    #print(path.kpts)

    bs = ph.band_structure(k)

    #print(bs)

    for i in range(len(bs)):
        for j in range(len(bs[i])):
            if(bs[i][j] > 0):
                Helm+= k_boltz*T*np.log(2*np.sinh(bs[i][j]/(2*k_boltz*T)))

    return(Helm/(N**3) + P*V + atoms.get_potential_energy()/N**3)

def toy_Helmholtz(atoms,calc,N,T,kvecs = (10, 10, 10)):

    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    ph.read()
    ph.clean()

    Helm = 0
    
    #path = atoms.cell.bandpath('GXULGK', npoints=100)
    k = monkhorst_pack(kvecs)
    #print(path.kpts)

    bs = ph.band_structure(k)

    for i in range(len(bs)):
        for j in range(len(bs[i])):
            if(bs[i][j] > 0):
                Helm+= k_boltz*T*np.log(2*np.sinh(bs[i][j]/(2*k_boltz*T)))
                
    return Helm/N**3 + atoms.get_potential_energy()/N**3


def main():
    a1: Atoms = bulk('Zr', 'hcp', a=3.24, c=5.18)
    a2: Atoms = bulk('Zr', 'bcc', a=3.24, cubic=True)
    calc = eamc.set_up_eam_calculator()
    a1.calc = calc
    a2.calc = calc
    free_energy_hcp = []
    free_energy_bcc = []
    lattice_constant_a = np.linspace(3.22, 3.26, 8)
    lattice_constant_c = np.linspace(5.175, 5.21, 8)
    count = 0
    for T in [50, 100, 200, 400, 600, 800, 1000, 1200]:
        a1.set_cell([lattice_constant_a[count], lattice_constant_a[count], lattice_constant_c[count], 90, 90, 120])
        a2.set_cell([lattice_constant_a[count], lattice_constant_a[count], lattice_constant_a[count]])
        free_energy_hcp.append(toy_gibbs(a1, calc, 6, 0.01, a1.get_volume(), T, kvecs=(12, 12, 8)))
        free_energy_bcc.append(toy_gibbs(a2, calc, 6, 0.01, a2.get_volume(), T, kvecs=(12, 12, 12)))
        #free_energy_hcp.append(toy_Helmholtz(a1, calc, 3, T))
        #free_energy_bcc.append(toy_Helmholtz(a2, calc, 3, T))
        count+=1
    plt.plot([50, 100, 200, 400, 600, 800, 1000, 1200], free_energy_bcc, label="bcc")
    plt.plot([50, 100, 200, 400, 600, 800, 1000, 1200], free_energy_hcp, label="hcp")
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    main()   
    
    
    



