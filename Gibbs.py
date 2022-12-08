from ase.build import bulk
from ase.calculators.emt import EMT
from ase.phonons import Phonons
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.dft.kpoints import monkhorst_pack
from ase.visualize import view
import ase.calculators.morse as morse
import eam_calculator as eamc
import meam_calculator as meamc
import sys

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
'''
def toy_gibbs(atoms,calc,N,P,V,T,kvecs = (6, 6, 6)):

    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    ph.read()
    ph.clean()

    Helm = 0
    entropy = 0
    #path = atoms.cell.bandpath('GXULGK', npoints=100)
    k = monkhorst_pack(kvecs)
    #print(path.kpts)

    bs = ph.band_structure(k)

    #print(bs)
    
    for i in range(len(bs)):
        #for j in range(len(bs[i])):
            #if(bs[i][j] > 0):
            Helm+= k_boltz*T*np.log(2*np.sinh(np.linalg.norm(bs[i], ord=1)/(2*k_boltz*T)))
            #entropy += -k_boltz*np.log(1-np.exp(-bs[i][j]/(k_boltz * T))) + bs[i][j]/((np.exp(-bs[i][j]/(k_boltz*T)-1)*T))
    return Helm/(N**3) + P*V + atoms.get_potential_energy()/N**3
'''

def toy_gibbs(atoms,calc,N,P,V,T,kvecs = (6, 6, 6)):

    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    ph.read(acoustic=True)
    ph.clean()

    k = monkhorst_pack(kvecs)
    bs = ph.band_structure(k, verbose=False)

    F1 = k_boltz*T*np.log(2*np.sinh(np.array(bs)/(2*k_boltz*T)))
    Helm = sum(map(sum, F1))
    return Helm/(N**3) + P*V + atoms.get_potential_energy()/N**3

def toy_Helmholtz(atoms,calc,N,T,kvecs = (10, 10, 10)):

    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    ph.read(acoustic=True)
    ph.clean()

    Helm = 0
    
    #path = atoms.cell.bandpath('GXULGK', npoints=100)
    k = monkhorst_pack(kvecs)
    #print(path.kpts)

    bs = ph.band_structure(k, verbose=False)

    for i in range(len(bs)):
        for j in range(len(bs[i])):
            if(bs[i][j] > 0):
                Helm+= k_boltz*T*np.log(2*np.sinh(bs[i][j]/(2*k_boltz*T)))
                
    return Helm/N**3 + atoms.get_potential_energy()/N**3


def main(argv):
    if argv[0] == "--meam":
        a1: Atoms = bulk('Zr', 'hcp', a=3.24, c=5.18)
        a2: Atoms = bulk('Zr', 'bcc', a=3.24, cubic=True)
        calc1 = meamc.setup_meam_calculator(a1)
        calc2 = meamc.setup_meam_calculator(a2)
        a1.calc = calc1
        a2.calc = calc2
        free_energy_hcp = []
        free_energy_bcc = []
        lattice_constant_a = np.linspace(3.22, 3.26, 9)
        lattice_constant_c = np.linspace(5.21, 5.28, 9)
        count = 0
        for T in [50, 100, 200, 400, 600, 800, 1000, 1200, 1400]:
            
            a1.set_cell([lattice_constant_a[count], lattice_constant_a[count], lattice_constant_c[count], 90, 90, 120])
            a2.set_cell([lattice_constant_a[count], lattice_constant_a[count], lattice_constant_a[count]])
            a1.set_positions(np.multiply(a1.get_positions(), np.array([lattice_constant_a[count] / 3.24, lattice_constant_a[count] / 3.24, lattice_constant_c[count]/5.18])))
            a2.set_positions(np.multiply(a2.get_positions(), np.array([lattice_constant_a[count] / 3.24, lattice_constant_a[count] / 3.24, lattice_constant_a[count]/3.24])))
            print(a1.get_positions())
            free_energy_hcp.append(toy_gibbs(a1, calc1, 8, 6.25e-3, a1.get_volume(), T, kvecs=(16, 16, 13)))
            free_energy_bcc.append(toy_gibbs(a2, calc2, 8, 6.25e-3, a2.get_volume(), T, kvecs=(16, 16, 16)))
            #free_energy_hcp.append(toy_Helmholtz(a1, calc, 3, T, kvecs=(6,6,4)))
            #free_energy_bcc.append(toy_Helmholtz(a2, calc, 3, T, kvecs=(6,6,6)))
            count+=1
        print(free_energy_bcc, free_energy_hcp)
        plt.plot([50, 100, 200, 400, 600, 800, 1000, 1200, 1400], free_energy_bcc, label="bcc")
        plt.plot([50, 100, 200, 400, 600, 800, 1000, 1200, 1400], free_energy_hcp, label="hcp")
        plt.legend()
        plt.xlabel("temperature K")
        plt.ylabel("Gibbs Energy per cell ev")
        plt.show()
    elif argv[0] == "--eam":
        a1: Atoms = bulk('Zr', 'hcp', a=3.24, c=5.18)
        a2: Atoms = bulk('Zr', 'bcc', a=3.24, cubic=True)
        calc1 = eamc.set_up_eam_calculator()
        a1.calc = calc1
        a2.calc = calc1
        free_energy_hcp = []
        free_energy_bcc = []
        lattice_constant_a = np.linspace(3.21, 3.26, 9)
        lattice_constant_c = np.linspace(5.16, 5.22, 9)
        count = 0
        for T in [50, 100, 200, 400, 600, 800, 1000, 1200, 1400]:
            
            a1.set_cell([lattice_constant_a[count], lattice_constant_a[count], lattice_constant_c[count], 90, 90, 120])
            a2.set_cell([lattice_constant_a[count], lattice_constant_a[count], lattice_constant_a[count]])
            #a1.set_positions(np.multiply(a1.get_positions(), np.array([lattice_constant_a[count] / 3.24, lattice_constant_a[count] / 3.24, lattice_constant_c[count]/5.18])))
            #a2.set_positions(np.multiply(a2.get_positions(), np.array([lattice_constant_a[count] / 3.24, lattice_constant_a[count] / 3.24, lattice_constant_a[count]/3.24])))
            print(a1.get_potential_energy())
            free_energy_hcp.append(toy_gibbs(a1, calc1, 8, 6.25e-3, a1.get_volume(), T, kvecs=(16, 16, 13)))
            free_energy_bcc.append(toy_gibbs(a2, calc1, 8, 6.25e-3, a2.get_volume(), T, kvecs=(16, 16, 16)))
            #free_energy_hcp.append(toy_Helmholtz(a1, calc, 3, T, kvecs=(6,6,4)))
            #free_energy_bcc.append(toy_Helmholtz(a2, calc, 3, T, kvecs=(6,6,6)))
            count+=1
        print(free_energy_bcc, free_energy_hcp)
        plt.plot([50, 100, 200, 400, 600, 800, 1000, 1200, 1400], free_energy_bcc, label="bcc")
        plt.plot([50, 100, 200, 400, 600, 800, 1000, 1200, 1400], free_energy_hcp, label="hcp")
        plt.legend()
        plt.xlabel("temperature K")
        plt.ylabel("Gibbs Energy per cell ev")
        plt.show()    
if __name__ == "__main__":
    main(sys.argv[1:])   
    
    
    



