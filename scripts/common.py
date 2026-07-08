"""
Common functions used in all plots and compute_implied_epsilon.py

"""

import numpy
from numpy import log, log10, sqrt, exp, pi, arange, logspace, linspace, r_, c_
# Don't judge me: 
from matplotlib.pyplot import *
from scipy.interpolate import UnivariateSpline as US
import hmf
from astropy.cosmology import FlatLambdaCDM

bo_color='#BF5700'
top_color='DodgerBlue'
my_cmap=get_cmap('plasma')

# Planck2020 parameters:
h0=0.6732
om0=0.3158
ob0h2=0.022383
ob0=ob0h2/h0**2
ns=0.96605
sig8=0.8120
t_cmb=2.7255
# f_bary=Omega_b/Omega_m
fbary=ob0h2/h0**2/om0


# set up cosmological model:
planck2020_model=FlatLambdaCDM(H0 = 100*h0, Om0=om0, Tcmb0 = t_cmb, Ob0 = ob0)

######################################################################
# Planck2020 + Early Dark Energy (EDE) parameters (Smith et al. 2022):
h0_ede=0.7483
om0_ede=0.287
ob0h2_ede=0.02278
ob0_ede=ob0h2_ede/h0_ede**2
ns_ede=1.003
sig8_ede=0.829/(om0_ede/0.3)**0.5
# f_bary=Omega_b/Omega_m
fbary_ede=ob0h2_ede/h0_ede**2/om0_ede

# set up cosmological model:
planck2020_ede=FlatLambdaCDM(H0 = 100*h0_ede, Om0=om0_ede,
                             Tcmb0 = t_cmb, Ob0 = ob0_ede)
######################################################################


######################################################################
# hmf parameters:
######################################################################
# standard min value to compute mass function (in log10):
# changed from 6 to 7
min_mval_log10=6
# min value to compute mass function (in log10) for nu calcs (should be low):
min_mval_nu_log10=-2
# max value to compute mass function (in log10) for nu calcs (should be hi):
max_mval_nu_log10=16.5

# comoving matter density (in h**2 Msun/Mpc^3)
rhom0=om0*277.8e9
# comoving matter density (in Msun/Mpc^3)
rhom0_noh=rhom0*h0**2
# comoving baryon density (in Msun/Mpc**3):  
rhob0_noh=rhom0_noh/om0*ob0
# comoving baryon density (in h**2 Msun/Mpc**3): 
rhob0=rhob0_noh/h0**2


######################################################################
# read in Labbe et al. data:
######################################################################
mydtype_mstar=numpy.dtype([('id', 'i4'), ('z', 'f4'),
                           ('z_low', 'f4'), ('z_hi', 'f4'),
                           ('z_xlow', 'f4'), ('z_xhi', 'f4'),  
                           ('log10_mstar', 'f4'),
                           ('log10_mstar_lo', 'f4'), ('log10_mstar_hi', 'f4'),
                           ('log10_mstar_xlo', 'f4'), ('log10_mstar_xhi', 'f4')])
# changed path to my data
mstar_data=numpy.loadtxt('/home/caltech-msanche3/data/jades_stellar_masses.dat',
                         dtype=mydtype_mstar, skiprows=1, # only need to skip header line #1
                         usecols=(0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

# sort by mass (descending)
print(mstar_data)
# dont need to sort data with one object
#my_argsort_mass=numpy.argsort(mstar_data['log10_mstar'])[::-1]
#mstar_data=mstar_data[my_argsort_mass]


mydtype_rhostar=numpy.dtype([('mstar', 'f4'),
                             ('rhostar', 'f4'), 
                             ('rhostar_lo', 'f4'), ('rhostar_hi', 'f4')])

# changed name to z_14 and only loaded one file since we are only focusing on one galaxy rn
# i also updated path to lead to my data
rhostar_data_z14=numpy.loadtxt('/home/caltech-msanche3/data/mass_cumul.dat', dtype=mydtype_rhostar)

# get redshifts of 2 most massive candidates
zz14=mstar_data['z']
print(zz14)
# SCRIPT LOOKS GOOD UP TO THIS POINT


######################################################################


# 4 pi steradians converted to arcmin^2:
arcmin2_fullsky=(180/pi)**2*4*pi*60**2

# JADES=GS area
jades_south_arcmin2=25.6096
# volume fraction is JADES GS area divided full sky volume
vol_frac=jades_south_arcmin2/arcmin2_fullsky
print(vol_frac)

# updated lows and highs to match the binning of our cum data
zlow_jades=13.0
# zmid_labbe=8.5 wont need this one
zhi_jades=15.0

# comoving volume between z=13 and z=15 in the JADES GS area in this cosmology;
# updated to z14 notation
vol_z14=(planck2020_model.comoving_volume(zhi_jades).value -
        planck2020_model.comoving_volume(zlow_jades).value)*vol_frac
print(vol_z14)

# same volumes in EDE:
# update to z14
vol_z14_ede=(planck2020_ede.comoving_volume(zhi_jades).value -
            planck2020_ede.comoving_volume(zlow_jades).value)*vol_frac


######################################################################
# get halo mass function data:
my_mf=hmf.MassFunction(Mmin=min_mval_log10, hmf_model='SMT',
                       cosmo_model=planck2020_model,
                       sigma_8=sig8, n=ns,
                       transfer_model=hmf.density_field.transfer_models.CAMB,
                       transfer_params={'extrapolate_with_eh':True})

# array of halo masses, for standard parameter inputs, this is a grid gong from min max log10 mass
mvals_noh=my_mf.m/h0


# a version that is appropriate for nu, where M_min (M_max) is lower (higher)
my_mf_nu=hmf.MassFunction(Mmin=min_mval_nu_log10, Mmax=max_mval_nu_log10,
                          hmf_model='SMT',
                          cosmo_model=planck2020_model,
                          sigma_8=sig8, n=ns,
                          transfer_model=hmf.density_field.transfer_models.CAMB,
                          transfer_params={'extrapolate_with_eh':True})
