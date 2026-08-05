import numpy as np
from flchem import flchem_evolve, flchem_init

# totaltime, timestep = 3.1556926e10, 3.1556926e10 #3.1556926e17, 7.8892315e16
# abundc, abundo = 1.4e-4,  3.2e-4
# dust_to_gas_ratio, dust_temp = 1e0, 10.
# G0, G0RP, cosmic_ray_ion_rate = 1.7, 0., 3e-17
# NH_ext, Z_atom = 1e20, 1.
# fshield_H2, fshield_CO, AV_mean, chi_mean = 1., 1., 0., 1.
# dl, divv = 0., 0.
# dens, cs = 2.1373888e-24, 1e5
# abh2, abhp, abco, abcp = .49, 1e-7, 1e-9, 0.

def calc_eint(dens,cs,gamma=1.6667):
    eint = dens*cs**2/(gamma*(gamma-1))
    return eint

# mH = 1.6726216370000000E-024
# abar = 1.4070615479999999

mH =   1.6726216370000000E-024
abar =   1.4074828305000000


totaltime, timestep = 3.1556926e10, 3.1556926e10 #3.1556926e17, 7.8892315e16
abundc, abundo = 1.4e-4,  3.2e-4
dust_to_gas_ratio, dust_temp = 1e0, 10.
G0, G0RP, cosmic_ray_ion_rate = 1.7, 0., 3e-17
NH_ext, Z_atom = 1e20, 1.
fshield_H2, fshield_CO, AV_mean, chi_mean = 1., 1., 0., 1.
dl, divv = 0., 0.
dens, cs = 2.1373888e-24, 1e5
abh2, abhp, abco, abcp = .49, 1e-7, 1e-9, 0.



#  time   315576000000.00000     
#  fshield_H2  0.90998533062480003     
#  fshield_CO  0.99336176001431920     
#  AV_mean  0.22000961256028617     
#  chi_mean  0.85791855966492092     
#  non_eq_abundances(ih2)   4.4370363056201866E-007
#  non_eq_abundances(ihp)   1.0240289394751549E-004
#  non_eq_abundances(ico)   8.4152660942134035E-018
#  dens   7.1070017170056703     
#  energy   6.6259635153971344E-012
#  ih2   4.7157830984661742E-007
#  ihp   1.1150618198888407E-004
#  ico   8.9457714208995610E-018
#  energy   6.3783181355411850E-012
#  temp   3939.1746275266246 

# time =  315576000000.00000     
# fshield_H2 =  0.90998533062480003     
# fshield_CO =  0.99336176001431920     
# AV_mean =  0.22000961256028617     
# chi_mean =  0.85791855966492092     
# abh2 =   4.4370363056201866E-007
# abhp =   1.0240289394751549E-004
# abco =   8.4152660942134035E-018
# dens = 7.1070017170056703 * mH * abar
# eint = 6.6259635153971344E-012

# time =  315576000000.00000     
# fshield_H2 =  0.90998533062480003     
# fshield_CO =  0.99336176001431920     
# AV_mean =  0.22000961256028617     
# chi_mean =  0.85791855966492092     
# abh2 =   4.4370363056201866E-007
# abhp =   1.0240289394751549E-004
# abco =   8.4152660942134035E-018
# dens = 7.1070017170056703 * mH * abar * (7.1070017170056703 / 7.1048744758159632)
# eint = 6.6259635153971344E-012

time = 315576000000.00000     
fshield_H2 =  0.90998532799172394     
fshield_CO =  0.99336175989629638     
AV_mean =  0.22000961256029070     
chi_mean =  0.85791855966491803     
abh2 =  4.4370364364480559E-007
abhp =  1.0240289394516484E-004
abco =  8.4152663425888911E-018
dens =  7.1048744517032354 * mH * abar 


time =  315576000000.00000     
fshield_H2 =  0.90278685696800465     
fshield_CO =  0.99310819859384181     
AV_mean =  0.21995985813416580     
chi_mean =  0.85795076371482704     
abh2 =  4.7285912078788135E-007
abhp =  1.0227777437299904E-004
abco =  8.9598445467324264E-018
dens =  7.1048744517032354 * mH * abar 
eint =  6.3551781809563165E-012



# 7.1070017170070097     * mH * abar #* (7.1070017170056703 / 7.1048744758159632)

print("dens before transform, before pychem call:", dens)
# dens = 1.6726197699478965E-023
# eint =  6.6259634867579533E-012



print("dens before pychem call:", dens)

# exit()
  



if __name__=='__main__':

    # eint = calc_eint(dens,cs)
    
    flchem_init(abundc, abundo, dust_to_gas_ratio, cosmic_ray_ion_rate, NH_ext, Z_atom)
    
    results = flchem_evolve(time, time,
                            G0, G0RP,
                            fshield_H2, fshield_CO, AV_mean, chi_mean,
                            dl, divv,
                            dens, eint,
                            abh2, abhp, abco, abcp)
                            
    print(results[:-2])
    print(results[-2])
    print(results[-1])

