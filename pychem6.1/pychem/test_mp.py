import numpy as np
from flchem import flchem_evolve, flchem_evolve_mp, flchem_init
from timeit import default_timer as timer
from multiprocessing import Pool, freeze_support

totaltime, timestep = 1e+10, 1e+10
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

    #dens_r = np.random.rand(10,8,8,8)*dens
    dens_r = np.ones((500,8,8,8))*dens
    eint_r = calc_eint(dens_r,cs)

    freeze_support()

    flchem_init(abundc, abundo, dust_to_gas_ratio, cosmic_ray_ion_rate, NH_ext, Z_atom)

    args = (totaltime, timestep, G0, G0RP,
             fshield_H2, fshield_CO, AV_mean, chi_mean, dl, divv, dens_r, eint_r,
             abh2, abhp, abco, abcp)
    
    start = timer()
    results = flchem_evolve_mp(*args)
    end = timer()

    cps = results[0].size/(end - start)
    print(f'Performance: {cps:.2f} cells/s')
