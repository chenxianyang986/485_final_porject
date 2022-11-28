from ase.build import bulk
from ase.calculators.emt import EMT
from ase.phonons import Phonons
import matplotlib.pyplot as plt
import numpy as np

h = 6.5821e-16
k = 8.617333262e-5


#METHOD 1: Discrete:
def gibbs(atoms,calc,N,P,V,kvecs = (10, 10, 10)):

    T = 10*i + 10

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


"""
def gibbs2(atoms,calc,N,P,V,kvecs = (10, 10, 10)):

    T = 10*i + 10

    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
    ph.run()

    ph.read()
    ph.clean()

    g_omega,dos = ph.dos(kvecs,1000,1e-3)

    #print(bs)

    funct = np.array(g_omega * 


    Helm = np.trapz(omega)

    return(Helm/(N**3) + P*V + atoms.get_potential_energy())

"""




    



