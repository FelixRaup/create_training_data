import numpy as np
from flchem import flchem_evolve, flchem_evolve_mp
from timeit import default_timer as timer

totaltime, timestep = 1.,1.#3.1556926e17, 7.8892315e16
abundc, abundo = 1.4e-4,  3.2e-4
dust_to_gas_ratio, dust_temp = 1e0, 10.
G0, cosmic_ray_ion_rate = 1.7, 3e-17
NH_ext, Z_atom = 1e20, 1e0
fshield_H2, fshield_CO, AV_mean, chi_mean = 1e0, 1e0, 0e0, 1e0
dl, divv = 0e0, 0e0
dens, cs = 2.1373888e-24, 1e5
abh2, abhp, abco, abcp = 4.9E-1, 1e-7, 1e-9, 0e0

def calc_eint(dens,cs,gamma=1.6667):
    eint = dens*cs**2/(gamma*(gamma-1))
    return eint

if __name__=='__main__':

    eint = calc_eint(dens,cs)
    
    
    
    start = timer()
    for i in range(100):
       results = flchem_evolve(totaltime, timestep,
                               abundc, abundo,
                               dust_to_gas_ratio, dust_temp,
                               G0, cosmic_ray_ion_rate,
                               NH_ext, Z_atom,
                               fshield_H2, fshield_CO, AV_mean, chi_mean,
                               dl, divv,
                               dens, eint,
                               abh2, abhp, abco, abcp)
    end = timer()
    print(results)
    print(end - start)
