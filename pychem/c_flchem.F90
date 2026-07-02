module flchem_interface
    use iso_c_binding, only: c_double, c_int
    implicit none    
    contains
    subroutine c_flinit(abundC, abundO &
                     ,dust_to_gas_ratio &
                     ,cosmic_ray_ion_rate &
                     ,NH_ext, Z_atom) bind(c) 

        real(c_double), intent(inout) :: abundC, abundO
        real(c_double), intent(inout) :: dust_to_gas_ratio
        real(c_double), intent(inout) :: cosmic_ray_ion_rate
        real(c_double), intent(inout) :: NH_ext, Z_atom

        call init_flchem(abundC, abundO &
                        ,dust_to_gas_ratio &
                        ,cosmic_ray_ion_rate &
                        ,NH_ext, Z_atom)
    end subroutine

    subroutine c_flchem(total_time, timestep &
                       ,G0, G0RP &
                       ,fshield_H2, fshield_CO, AV_mean, chi_mean &
                       ,dl, divv &
                       ,density, internal_energy &
                       ,abh2, abhp, abco, abcp &
                       ,temperature, dust_temp, tend &
                       ,istatus &
                       ,cool_rates, cool_rates_chem) bind(c) 

        real(c_double), intent(inout) :: total_time, timestep
        real(c_double), intent(inout) :: G0, G0RP
        real(c_double), intent(inout) :: fshield_H2, fshield_CO, AV_mean, chi_mean
        real(c_double), intent(inout) :: dl, divv
        real(c_double), intent(inout) :: density, internal_energy
        real(c_double), intent(inout) :: abh2, abhp, abco, abcp
        !real(c_double), dimension(??), intent(inout) :: pilc, pelc ! Hack: Not using those for now
        real(c_double), intent(inout) :: temperature, dust_temp, tend
        integer(c_int), intent(inout) :: istatus
        real(c_double), dimension(30), intent(inout) :: cool_rates
        real(c_double), dimension(8), intent(inout) :: cool_rates_chem

        call test_flchem(total_time, timestep &
                        ,G0, G0RP &
                        ,fshield_H2, fshield_CO, AV_mean, chi_mean &
                        ,dl, divv &
                        ,density, internal_energy &
                        ,abh2, abhp, abco, abcp &
#ifdef TR_ONTHESPOT
                        ,pilc, pelc &
#endif
                        ,temperature, dust_temp, tend &
                        ,istatus &
                        ,cool_rates, cool_rates_chem)

    end subroutine
end module
