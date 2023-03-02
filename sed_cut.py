
import os
import shutil
from astropy.io import ascii, fits
import numpy as np

object_file = 'fullcat307_315.fits'
objects_path = os.path.abspath('input_test/' + object_file)
output_folder = os.path.join(os.path.abspath(''), 'output_test')
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)
output_file = os.path.abspath(f'cut_{object_file}')
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
pm_dec_cond = abs(pm_dec_R)<2
data_pm_cut = data_pmra_cut[~pm_dec_cond]
print('----------TABLE ROWS AFTER pm_dec CUTTING----------')
print(len(data_pm_cut))

# ----------------------CONDITION 3---------------------------
mag_g = data_pm_cut['mag_auto_g']
mag_g_cond1 = mag_g==0
mag_g_cond2 = mag_g==99
data_g_cut = data_pm_cut[mag_g_cond1 | mag_g_cond2 ]
# data_g_delve_cut = data_g_delve_cut1 & data_g_delve_cut2
print('----------TABLE ROWS AFTER mag_auto_g CUTTING----------')
print(len(data_g_cut))

# ----------------------CONDITION 4 optional---------------------------
wise1 = data_g_cut['w1mpro']
wise2 = data_g_cut['w2mpro']
#wise2_cond = wise2!=0
wise_cond1 = (wise1+2.699)-(wise2+3.399)<= 0.6
wise_cond2 = (wise1+2.699)-(wise2+3.399)>= -0.6
data_wise_cut = data_g_cut[wise_cond1 & wise_cond2]
print('----------TABLE ROWS AFTER wise CUTTING----------')
print(len(data_wise_cut))

fits.writeto(output_file, data_wise_cut, overwrite=True)