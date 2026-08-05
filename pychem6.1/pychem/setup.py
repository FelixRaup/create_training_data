from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

flash_flags = '-DCHEMISTRYNETWORK=5 -g -fbacktrace -fcheck=all -std=legacy -O2 -fdefault-real-8 -ffree-line-length-0 -fallow-argument-mismatch -fdefault-double-8 -fbounds-check'

flash77_modules = ('coolinmo',
                 'cheminmo',
                 'spline',
                 'cool_func',
                 'photoinit_ism',
                 'dvode',
                 'evolve_abundances',
                 'rate_eq',
                 'jac',
                 'cool_util',
                 'const_rates',
                 'validate_rates',
                 'test_flchem',
                 'init_flchem',
                 'tr_rpGetKappa',
                 )
fl77_files = [f'{mn}.F' for mn in flash77_modules]
fl77_objects = [f'{mn}.o' for mn in flash77_modules]

#flash90_modules = ('tr_rpGetKappa',
#                  )
#fl90_files = [f'{mn}.F90' for mn in flash90_modules]
#fl90_objects = [f'{mn}.o' for mn in flash90_modules]

ifc_files = ['c_flchem.F90',]
ifc_objects = ['c_flchem.o',]

# compile the Flash Fortran 77 modules without linking
for fn,fo in zip(fl77_files,fl77_objects):
      fortran_mod_comp = f'gfortran {fn} -c -o {fo} -fPIC {flash_flags}'
#      print(fortran_mod_comp)
      system(fortran_mod_comp)
# compile the Flash Fortran 90 modules without linking
#for fn,fo in zip(fl90_files,fl90_objects):
#      fortran_mod_comp = f'gfortran {fn} -c -o {fo} -fPIC {flash_flags}'
#      print(fortran_mod_comp)
#      system(fortran_mod_comp)
# compile the Fortran interface module without linking
for fn,fo in zip(ifc_files,ifc_objects):
      comline = f'gfortran {fn} -c -o {fo} -O3 -fPIC'
#      print(comline)
      system(comline)

ext_modules = [Extension('flchem',# module name
                         ['flchem.pyx'], # source file
                         #language="c++",
                         extra_compile_args=['-fPIC', '-O3', '-fopenmp',],
                         libraries = ['gfortran', 'm', 'mvec'],
                         extra_link_args=fl77_objects+ifc_objects+['-fopenmp','-lgfortran',"-Wl,--no-as-needed", "-lmvec"]
                         )]

setup(name = 'pychem',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
