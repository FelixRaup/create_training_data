import numpy as np
from flchem import flchem_evolve, flchem_init

totaltime, timestep = 1., 1. #3.1556926e17, 7.8892315e16
abundc, abundo = 1.4e-4,  3.2e-4
dust_to_gas_ratio, dust_temp = 1e0, 10.
G0, G0RP, cosmic_ray_ion_rate = 1.7, 0., 3e-17
NH_ext, Z_atom = 1e20, 1.
fshield_H2, fshield_CO, AV_mean, chi_mean = 1., 1., 0., 1.
dl, divv = 0., 0.
dens, cs = 2.1373888e-24, 1e5
abh2, abhp, abco, abcp = .49, 1e-7, 1e-9, 0.

def calc_eint(dens,cs,gamma=1.6667):
    eint = dens*cs**2/(gamma*(gamma-1))
    return eint

if __name__=='__main__':

    eint = calc_eint(dens,cs)
    
    flchem_init(abundc, abundo, dust_to_gas_ratio, cosmic_ray_ion_rate, NH_ext, Z_atom)
    
    results = flchem_evolve(totaltime, timestep,
                            G0, G0RP,
                            fshield_H2, fshield_CO, AV_mean, chi_mean,
                            dl, divv,
                            dens, eint,
                            abh2, abhp, abco, abcp)
                            
    print(results[:-2])
    print(results[-2]/dens)
    print(results[-1]/dens)

