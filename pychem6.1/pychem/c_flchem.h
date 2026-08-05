extern void c_flchem(double* total_time, double* timestep,
                     double* G0, double* G0RP,
                     double* fshield_H2, double* fshield_CO, double* AV_mean, double* chi_mean,
                     double* dl, double* divv,
                     double* density, double* internal_energy,
                     double* abh2, double* abhp, double* abco, double* abcp,
                     double* temperature, double* dust_temp, double* tend,
                     long* istatus,
                     double* cool_rates, double* cool_rates_chem);