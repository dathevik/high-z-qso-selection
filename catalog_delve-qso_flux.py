#!/usr/bin/python

import numpy as np
import os
from astropy.io import fits
from astropy.io import ascii
import astropy.units as u # to define units

############################################################################################################
############################################################################################################

dir = os.path.abspath('input_test/')
############################################################################################################
############################################################################################################

object_file = '/cut_fullcat345_360.fits'
file = fits.open(dir+object_file)  # open a FITS file

hdr_o = file[0].header
#list(hdr.keys())

t = file[1].data
cols = t.columns
#cols.names
#print(t.columns)
#print(cols.names)


for col in cols.names[1:21]:
	if hasattr(t[col], "mask"):
		t[col][t[col].mask] = np.nan
	w = np.where((t[col] == -99.0) | (t[col] == 9999.0) | (t[col] < 0.) | (t[col] > 50.))
	if len(w[0]) > 0:
		t[col][w] = np.nan


file.writeto(dir+f"{object_file}_nan.fits", overwrite=True)  ## create the new mag nan fits table
file.close()


##############

f = fits.open(dir+f'{object_file}_nan.fits')  # open a FITS file

hdr = f[0].header
#list(hdr.keys())

tbdata = f[1].data
cols = tbdata.columns
cols.names


id_list = tbdata['quick_object_id'][:]
print('ids: ')
print(len(id_list))

ra = tbdata['ra'][:]
dec = tbdata['dec'][:]

#### reddhift photometrico
z_1 = np.ones((len(id_list)))
z = z_1*-1.
print (z)
print (len(z))



#### From AB mag to  flux [Jy]: AB=2.5*(23-log10(F/Jy))-48.6


def mag_to_flux(mag, err):	
	print(mag, err)
	flux_mJy = (mag*u.ABmag).to('mJy')

	mag_2 = mag + err
	print(mag_2)
	flux_err_mJy = (mag_2*u.ABmag).to('mJy')
	err_final = np.abs(flux_mJy - flux_err_mJy)
	
	return flux_mJy, err_final


####################################  SDSS  ##########################################
###### u
#uflux_mJy_sdss, uerr_final_sdss = mag_to_flux(tbdata['petromag_u_sdss'][:], tbdata['petromagerr_u_sdss'][:])
###### g
#gflux_mJy_sdss, gerr_final_sdss = mag_to_flux(tbdata['petromag_g_sdss'][:], tbdata['petromagerr_g_sdss'][:])
###### r
#rflux_mJy_sdss, rerr_final_sdss = mag_to_flux(tbdata['petromag_r_sdss'][:], tbdata['petromagerr_r_sdss'][:])
###### i
#iflux_mJy_sdss, ierr_final_sdss = mag_to_flux(tbdata['petromag_i_sdss'][:], tbdata['petromagerr_i_sdss'][:])
###### z
#zflux_mJy_sdss, zerr_final_sdss = mag_to_flux(tbdata['petromag_z_sdss'][:], tbdata['petromagerr_z_sdss'][:])


####################################  DES  ##########################################
###### g
#gflux_mJy_des, gerr_final_des = mag_to_flux(tbdata['petromag_g_des'][:], tbdata['petromagerr_g_des'][:])
###### r
#rflux_mJy_des, rerr_final_des = mag_to_flux(tbdata['petromag_r_des'][:], tbdata['petromagerr_r_des'][:])
###### i
#iflux_mJy_des, ierr_final_des = mag_to_flux(tbdata['petromag_i_des'][:], tbdata['petromagerr_i_des'][:])
###### z
#zflux_mJy_des, zerr_final_des = mag_to_flux(tbdata['petromag_z_des'][:], tbdata['petromagerr_z_des'][:])


####################################  DELVE  ##########################################
###### g
gflux_mJy_delve, gerr_final_delve = mag_to_flux(tbdata['mag_auto_g'][:], tbdata['magerr_auto_g'][:])
###### r
rflux_mJy_delve, rerr_final_delve = mag_to_flux(tbdata['mag_auto_r'][:], tbdata['magerr_auto_r'][:])
###### i
iflux_mJy_delve, ierr_final_delve = mag_to_flux(tbdata['mag_auto_i'][:], tbdata['magerr_auto_i'][:])
###### z
zflux_mJy_delve, zerr_final_delve = mag_to_flux(tbdata['mag_auto_z'][:], tbdata['magerr_auto_z'][:])


####################################  VHS  ##########################################

###### Y mag
Y_ABmag = tbdata['ypetromag'][:] + 0.6
yflux_mJy_vhs, yerr_final_vhs = mag_to_flux(Y_ABmag, tbdata['ypetromagerr'][:])

###### J mag
J_ABmag = tbdata['jpetromag'][:] + 0.916
jflux_mJy_vhs, jerr_final_vhs = mag_to_flux(J_ABmag, tbdata['jpetromagerr'][:])

###### H mag
H_ABmag = tbdata['hpetromag'][:] + 1.366
hflux_mJy_vhs, herr_final_vhs = mag_to_flux(H_ABmag, tbdata['hpetromagerr'][:])

###### Ks mag
KS_ABmag = tbdata['kspetromag'][:] + 1.827
ksflux_mJy_vhs, kserr_final_vhs = mag_to_flux(KS_ABmag, tbdata['kspetromagerr'][:])


####################################  UKIDSS  ##########################################
###### Y mag
#yflux_mJy_ukidss, yerr_final_ukidss = mag_to_flux(tbdata['petromag_Y_ukidss'][:], tbdata['petromagerr_Y_ukidss'][:])
###### J mag
#jflux_mJy_ukidss, jerr_final_ukidss = mag_to_flux(tbdata['petromag_J_ukidss'][:], tbdata['petromagerr_J_ukidss'][:])
###### H mag
#hflux_mJy_ukidss, herr_final_ukidss = mag_to_flux(tbdata['petromag_H_ukidss'][:], tbdata['petromagerr_H_ukidss'][:])
###### Ks mag
#kflux_mJy_ukidss, kerr_final_ukidss = mag_to_flux(tbdata['petromag_K_ukidss'][:], tbdata['petromagerr_K_ukidss'][:])


####################################  WISE  ##########################################
###### W1
W1_ABmag = tbdata['w1mpro'][:] + 2.699
w1flux_mJy, w1err_final = mag_to_flux(W1_ABmag, tbdata['w1sigmpro'][:])

###### W2
W2_ABmag = tbdata['w2mpro'][:] + 3.339
w2flux_mJy, w2err_final = mag_to_flux(W2_ABmag, tbdata['w2sigmpro'][:])
###### W3
#w3flux_mJy, w3err_final = mag_to_flux(tbdata['intmag_w3_allwise'][:], tbdata['intmagerr_w3_catwise'][:])
###### W4
#w4flux_mJy, w4err_final = mag_to_flux(tbdata['intmag_w4_allwise'][:], tbdata['intmagerr_w4_catwise'][:])

############# W3 upper limits ########
 
#new_W3, new_W3err = [],[] # w3flux_mJy, w3err_final
#for i in range(len(w3flux_mJy)):

#	if (w3flux_mJy[i] > 0. and w3err_final[i] > 0.): 
#		print('detection:   W3 = '+str(w3flux_mJy[i])+ ', W3_err = '+str(w3err_final[i]))
#		new_W3.append(float(w3flux_mJy[i].value))
#		new_W3err.append(float(w3err_final[i].value))
		
#	elif (w3flux_mJy[i] > 0. and str(w3err_final[i].value) == 'nan'): 
#		print('upper limit:   W3 = '+str(w3flux_mJy[i])+ ', W3_err = '+str(w3err_final[i]))
#		new_flux_err = w3flux_mJy[i].value*-1.
#		new_W3.append(float(w3flux_mJy[i].value))
#		new_W3err.append(float(new_flux_err))	
		
#	elif str(w3flux_mJy[i].value) == 'nan':
#		print('no information:   W3 = '+str(w3flux_mJy[i])+ ', W3_err = '+str(w3err_final[i]))
#		new_W3.append(float(w3flux_mJy[i].value))
#		new_W3err.append(float(w3err_final[i].value))
		
#	else: continue
		
########### W4 upper limits #########
 
#new_W4, new_W4err = [],[]
#for j in range(len(w4flux_mJy)):

#	if (w4flux_mJy[j] > 0. and w4err_final[j] > 0.): 
#		print('detection:   W4 = '+str(w4flux_mJy[j])+ ', W4_err = '+str(w4err_final[j]))
#		new_W4.append(float(w4flux_mJy[j].value))
#		new_W4err.append(float(w4err_final[j].value))
		
#	elif (w4flux_mJy[j] > 0. and str(w4err_final[j].value) == 'nan'): 
#		print('upper limit:   W4 = '+str(w4flux_mJy[j])+ ', W4_err = '+str(w4err_final[j]))
#		new_flux_err = w4flux_mJy[j].value*-1.
#		new_W4.append(float(w4flux_mJy[j].value))
#		new_W4err.append(float(new_flux_err))	
		
#	elif str(w4flux_mJy[j].value) == 'nan':
#		print('no information:   W4 = '+str(w4flux_mJy[j])+ ', W4_err = '+str(w4err_final[j]))
#		new_W4.append(float(w4flux_mJy[j].value))
#		new_W4err.append(float(w4err_final[j].value))
		
#	else: continue



####################################  GALEX  ##########################################
###### NUV
#nuvflux_mJy, nuverr_final = mag_to_flux(tbdata['kronmag_NUV_galex'][:], tbdata['kronmagerr_NUV_galex'][:])
###### FUV
#fuvflux_mJy, fuverr_final = mag_to_flux(tbdata['kronmag_FUV_galex'][:], tbdata['kronmagerr_FUV_galex'][:])



#### output file
#uflux_mJy_sdss, uerr_final_sdss, gflux_mJy_sdss, gerr_final_sdss, rflux_mJy_sdss, rerr_final_sdss, iflux_mJy_sdss, ierr_final_sdss, zflux_mJy_sdss, zerr_final_sdss, 
#gflux_mJy_des, gerr_final_des, rflux_mJy_des, rerr_final_des, iflux_mJy_des, ierr_final_des, zflux_mJy_des, zerr_final_des, 
			
			
#'u_prime_sdss', 'u_prime_err_sdss', 'g_prime_sdss', 'g_prime_err_sdss', 'r_prime_sdss', 'r_prime_err_sdss', 'i_prime_sdss', 'i_prime_err_sdss', 'z_prime_sdss', 'z_prime_err_sdss',
#'g_prime_des', 'g_prime_err_des', 'r_prime_des', 'r_prime_err_des', 'i_prime_des', 'i_prime_err_des', 'z_prime_des', 'z_prime_err_des',
            	   

ascii.write([id_list, ra, dec, z, 
			gflux_mJy_delve, gerr_final_delve, rflux_mJy_delve, rerr_final_delve, iflux_mJy_delve, ierr_final_delve, zflux_mJy_delve, zerr_final_delve, 
			yflux_mJy_vhs, yerr_final_vhs, jflux_mJy_vhs, jerr_final_vhs, hflux_mJy_vhs, herr_final_vhs, ksflux_mJy_vhs, kserr_final_vhs,
			w1flux_mJy, w1err_final, w2flux_mJy, w2err_final], 
            dir+f'{object_file}_fluxes.dat',
            names=['#id', 'ra', 'dec', 'redshift', 
            	   'g_prime_delve', 'g_prime_err_delve', 'r_prime_delve', 'r_prime_err_delve', 'i_prime_delve', 'i_prime_err_delve', 'z_prime_delve', 'z_prime_err_delve',
            	   'vista.vircam.Y_vhs', 'vista.vircam.Y_err_vhs', 'vista.vircam.J_vhs', 'vista.vircam.J_err_vhs', 'vista.vircam.H_vhs', 'vista.vircam.H_err_vhs', 'vista.vircam.Ks_vhs', 'vista.vircam.Ks_err_vhs',  
                   'WISE1', 'WISE1_err', 'WISE2', 'WISE2_err'],  
                   overwrite=True)
         
         # yflux_mJy_ukidss, yerr_final_ukidss, jflux_mJy_ukidss, jerr_final_ukidss, hflux_mJy_ukidss, herr_final_ukidss, kflux_mJy_ukidss, kerr_final_ukidss,
         # new_W3, new_W3err, new_W4, new_W4err,            
		 # nuvflux_mJy, nuverr_final, fuvflux_mJy, fuverr_final], 
		 
		 # 'ukirt.Y', 'ukirt.Y_err', 'ukirt.J', 'ukirt.J_err', 'ukirt.H', 'ukirt.H_err', 'ukirt.Ks', 'ukirt.Ks_err',
		 # 'WISE3', 'WISE3_err', 'WISE4', 'WISE4_err',
		 # 'galex.NUV', 'galex.NUV_err', 'galex.FUV', 'galex.FUV_err'
		 
		 
		 
		 