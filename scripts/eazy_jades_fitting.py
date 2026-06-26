# code made specifically for HellBender cluster fitting
#import pandas as pd
import eazy
import os
import eazy.photoz
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from eazy import filters, utils, templates

# path to JADES catalog on HellBender 
cat_name = '/home/caltech-msanche3/data/JADES_fits_data.fits'

# path to translate file on HellBender
trans_file = '/home/caltech-msanche3/data'

jades_table = Table.read(cat_name, hdu='KRON_CONV',format='fits')
params = {}
print("Loading parameters...")
# 1. Tell it where your data is
params['CATALOG_FILE'] = jades_table
#params['CATALOG_FORMAT'] = 'fits'
params['EXT_NUMBER'] = 9
params['SYS_ERR'] = 0.03

# 2. Tell it what templates to use to fit the light
params['TEMPLATES_FILE'] = "/home/caltech-msanche3/data/miniconda/envs/eazyenv2/lib/python3.13/site-packages/eazy/data/eazy-photoz/templates/spline_templates_v2/tweak_spline.param"

# 3. Choose a base name for the output files it saves
# path to jades data output directory
# ^^ let's not do this for now and just name the file, we can see where it is directed to once code runs
params['MAIN_OUTPUT_FILE'] = 'jades_run'

print("Loaded parameters!")
print('Initializing fitting program...')
# initialize the fitting program
#add path to translate file
trans_file = '/home/caltech-msanche3/data/jades_trans_file.csv'
ez = eazy.photoz.PhotoZ(
    param_file=None, # don't call default file
    translate_file=trans_file, # call translate file
    zeropoint_file=None, # no zero points file for now?
    params=params, # call your own params, eazy will use default for those not defined
    load_prior=False, # uh
    load_products=False # uh
    
)

print('Fitting program initialized successfully!')

ez.set_sys_err(positive=True) # should I change this to False?

# 2. Fit the ENTIRE catalog
print('Fitting catalog in chunks...')
# total objects: 304366
# split into chunks of 10000
for i in range(0, 304366, 10000):
    print(f'Fitting chunk {i//10000 + 1}...')
    chunk = ez.idx[i : i + 10000]
    ez.fit_catalog(chunk,n_proc=16)
print('Catalog fitted successfully!')

print('Generating standard output...')
# Derived parameters (z params, RF colors, masses, SFR, etc.)
zout, hdu = ez.standard_output(
    simple=False, 
    rf_pad_width=0.5,
    rf_max_err=2, 
    prior=True,
    beta_prior=True, 
    absmag_filters=[], 
    extra_rf_filters=[]
)
# 'zout' also saved to [MAIN_OUTPUT_FILE].zout.fits

print('Standard output generated successfully!')
print('All steps completed successfully!')
print('Output saved to [MAIN_OUTPUT_FILE].zout.fits')
