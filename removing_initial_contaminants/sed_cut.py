import os
import shutil
import numpy as np
from astropy.io import fits

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Input and output paths relative to script directory
object_file = 'fullcat_viking_second_area.fits'
objects_path = os.path.join(script_dir, object_file)
# Get the directory of the input file
input_dir = os.path.dirname(objects_path)
# Create output file in the same directory as input
output_file = os.path.join(input_dir, f'cut_{object_file}')

print("READING OBJECTS:", objects_path, "FILE")
data_obj = fits.open(objects_path)
data_objects = data_obj[1].data
cols = data_objects.columns
print(data_objects)
ext_z_class = data_objects['extended_class_z']

print('----------TABLE ROWS BEFORE CUTTING----------')
print(len(data_objects))

# ----------------------CONDITION 1---------------------------
z_cond = ext_z_class==3
data_z_cut = data_objects[~z_cond]
print('----------TABLE ROWS AFTER ext_z CUTTING----------')
print(len(data_z_cut))

# ----------------------CONDITION 2---------------------------
pm_ra_err = data_z_cut['pmra_error']
pm_ra = data_z_cut['pmra']
pm_ra_R = pm_ra / pm_ra_err
pm_ra_cond=abs(pm_ra_R)>2
data_pmra_cut = data_z_cut[~pm_ra_cond]
print('----------TABLE ROWS AFTER pm_ra CUTTING----------')
print(len(data_pmra_cut))
pm_dec = data_pmra_cut['pmdec']
pm_dec_err = data_pmra_cut['pmdec_error']
pm_dec_R = pm_dec/pm_dec_err
pm_dec_cond = abs(pm_dec_R)>2
data_pm_cut = data_pmra_cut[~pm_dec_cond]
print('----------TABLE ROWS AFTER pm_dec CUTTING----------')
print(len(data_pm_cut))

# ----------------------CONDITION 3---------------------------
mag_g = data_pm_cut['mag_auto_g']
mag_g_cond1 = mag_g==0
mag_g_cond2 = mag_g==99
data_g_cut = data_pm_cut[mag_g_cond1 | mag_g_cond2]
print('----------TABLE ROWS AFTER mag_auto_g CUTTING----------')
print(len(data_g_cut))

# ----------------------CONDITION 4 optional---------------------------
wise1 = data_g_cut['w1mpro']
wise2 = data_g_cut['w2mpro']
wise1_cond0 = np.isnan(wise1)
wise2_cond0 = np.isnan(wise2)
wise_cond1 = (wise1+2.699) - (wise2+3.399) <= 0.6
wise_cond2 = (wise1+2.699) - (wise2+3.399) >= -0.6

data_wise_cut = data_g_cut[wise_cond1 & wise_cond2 | wise1_cond0 | wise2_cond0]
print('----------TABLE ROWS AFTER wise CUTTING----------')
print(len(data_wise_cut))

fits.writeto(output_file, data_wise_cut, overwrite=True) 