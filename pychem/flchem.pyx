cimport numpy as np
import numpy as np
from numpy cimport ndarray as ar
from multiprocessing import Pool, freeze_support
cimport cython
from cython.parallel cimport parallel, prange

cdef extern from "c_flchem.h":
    void c_flchem(double* simtime, double* simdt,
                  double* G0, double* G0RP,
                  double* fshield_H2, double* fshield_CO, double* AV_mean, double* chi_mean,
                  double* dl, double* divv,
                  double* dens, double* eint,
                  double* abh2, double* abhp, double* abco, double* abcp,
                  double* temperature, double* dust_temp, double* tend,
                  long* istatus,
                  double* cool_rates, double* cool_rates_chem) noexcept nogil

cdef extern from "c_flinit.h":
    void c_flinit(double* abundC, double* abundO,
                  double* dust_to_gas_ratio,
                  double* cosmic_ray_ion_rate,
                  double* NH_ext, double* Z_atom) noexcept nogil


def outershape(*ar_):
    dshapes = [np.shape(a) for a in ar_]
    dndim = [np.ndim(a) for a in ar_]
    ndim_outer = np.max(dndim)
    oshape = list()
    for i in range(ndim_outer):
        s_i = [ds[::-1][i] for ds in dshapes if len(ds)>i]
        smax_i = np.max(s_i)
        oshape.append(smax_i)
    return tuple(oshape[::-1])

cpdef project(a,dshape):
    arr = np.array(a)
    ones = np.ones(dshape, dtype=arr.dtype)
    ap = arr*ones
    return ap

def argzip(*ar_):
    oshape = outershape(*ar_)
    proj_linshape = lambda d: project(d,oshape).ravel()
    ar_oshape = map(proj_linshape,ar_)
    ar_linshape = zip(*ar_oshape)
    return tuple(ar_linshape)
    
def flchem_evolve_mp(simtime, simdt,
                  G0, G0RP,
                  fshield_H2, fshield_CO, AV_mean, chi_mean,
                  dl, divv,
                  dens, eint,
                  abh2, abhp, abco, abcp,
                  ):
                  
    args = (simtime, simdt,
                  G0, G0RP,
                  fshield_H2, fshield_CO, AV_mean, chi_mean,
                  dl, divv,
                  dens, eint,
                  abh2, abhp, abco, abcp)

    # def flchem_evolve_mp(*args) would also work, but above argument list
    # was left in in order to check and document the expected arguments

    freeze_support()

    oshape = outershape(*args)
    linargs = argzip(*args)

    with Pool() as mp:
        result = mp.starmap(flchem_evolve, linargs)

    result_arrays = [np.array(d) for d in zip(*result)]

    scalars = [a.reshape(oshape) for a in result_arrays[:-2]]
    cool_rates = result_arrays[-2].T.reshape((30,)+oshape)
    cool_rates_chem = result_arrays[-1].T.reshape((8,)+oshape)

    return *result_arrays[:-2], cool_rates, cool_rates_chem


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef flchem_init(double abundC, double abundO,
                  double dust_to_gas_ratio,
                  double cosmic_ray_ion_rate,
                  double NH_ext, double Z_atom,
                  ) noexcept:

    c_flinit(&abundC, &abundO,
            &dust_to_gas_ratio,
            &cosmic_ray_ion_rate,
            &NH_ext, &Z_atom)

    return None



@cython.wraparound(False)
@cython.boundscheck(False)
cpdef flchem_evolve(double simtime, double simdt,
                  double G0, double G0RP,
                  double fshield_H2, double fshield_CO, double AV_mean, double chi_mean,
                  double dl, double divv,
                  double dens, double eint,
                  double abh2, double abhp, double abco, double abcp,
                  ) noexcept:


    cdef ar[double, ndim=1] cool_rates = np.empty(30, order='F')
    cdef ar[double, ndim=1] cool_rates_chem = np.empty(8, order='F')

    dust_temp = 10.
    temperature = 0.
    tend = 0.
    istatus = 4711

    c_flchem(&simtime, &simdt,
               &G0, &G0RP,
               &fshield_H2, &fshield_CO, &AV_mean, &chi_mean,
               &dl, &divv,
               &dens, &eint,
               &abh2, &abhp, &abco, &abcp,
               &temperature, &dust_temp, &tend,
               &istatus,
               <double*> cool_rates.data, <double*> cool_rates_chem.data)

    result = (    tend, # istatus,
                  abh2, abhp, abco, abcp,
                  temperature, dust_temp, eint, istatus,
                  cool_rates, cool_rates_chem)
    
    return result
