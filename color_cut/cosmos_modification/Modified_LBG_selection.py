
""""This activity is intended to correct intergalactic medium (IGM) absorption for our selection criteria.
 We are using the Lyman break galaxy (LBG) composite spectrum from Shaply et al. (2003); 
 in their study, they used the IGM absorption for redshift 3, 
 to apply the Madau (1995) Lyman alpha forest at targeted redshift, 
 we removed the IGM absorption from the LBG composite Spectrum at z~3 and we applied 
the IGM absorption at the Hot DOG redshift, by taking into account the IGM absorption from Madau (1995)."""


import numpy as np
import matplotlib.pyplot as plt
from synphot import SourceSpectrum, Observation, units
from synphot.spectrum import SpectralElement
from synphot.models import Empirical1D
from synphot import SourceSpectrum, etau_madau

#Load the spectrum.
#data = np.loadtxt("qsogen_sed_emlines0.0_ebv+000_z+500.dat")
data = np.loadtxt("qsogen_sed_emlines0.0_ebv+000_z+500.dat")
LBG = SourceSpectrum(Empirical1D, points=data[:,0], lookup_table=data[:,1]*units.FNU, keep_neg=True)


#HSC_g = np.loadtxt("Subaru_HSC_g.dat")
HSC_r = np.loadtxt("Subaru_HSC.r.dat")
HSC_i = np.loadtxt("Subaru_HSC.i.dat")
HSC_z = np.loadtxt("Subaru_HSC.z.dat")


#imacs_g = np.loadtxt("interpolated_sdss_filter_g.dat")
imacs_r = np.loadtxt("CTIO_DECam.r.dat")
imacs_i = np.loadtxt("CTIO_DECam.i.dat")
imacs_z = np.loadtxt("CTIO_DECam.z.dat")

#3.427<z< 3.827
#3.712<z< 4.112

z=5

# LBG.z =z
# wave_shi = data[:,0] * (1 + z)
# print(wave_shi)
# extcurve = etau_madau(wave_shi, z)
# tau_eff = -np.log(extcurve(data[:,0]*(1+z)))
# kappa = 1
# tau_use = kappa*tau_eff
# lbg_fnu_with_igm = data[:,1] * np.exp(-tau_use) * units.FNU
# LBG_with_IGM = SourceSpectrum(Empirical1D, points=wave_shi, lookup_table=lbg_fnu_with_igm, keep_neg=True)



# bp_HSC_g = SpectralElement(Empirical1D, points=HSC_g[:,0], lookup_table=HSC_g[:,1])
bp_HSC_r = SpectralElement(Empirical1D, points=HSC_r[:,0], lookup_table=HSC_r[:,1])
bp_HSC_i = SpectralElement(Empirical1D, points=HSC_i[:,0], lookup_table=HSC_i[:,1])
bp_HSC_z = SpectralElement(Empirical1D, points=HSC_z[:,0], lookup_table=HSC_z[:,1])

# obs_hg = Observation(LBG_with_IGM, bp_HSC_g, force = 'extrap')#, binset='array-like')
obs_hr = Observation(LBG, bp_HSC_r, force= 'extrap')#, binset='array-like')
obs_hi = Observation(LBG, bp_HSC_i, force = 'extrap')#,  binset='array-like')
obs_hz = Observation(LBG, bp_HSC_z, force = 'extrap')#,  binset='array-like')


#bp_imacs_g = SpectralElement(Empirical1D, points=imacs_g[:,0], lookup_table=imacs_g[:,1])
bp_imacs_r = SpectralElement(Empirical1D, points=imacs_r[:,0], lookup_table=imacs_r[:,1])
bp_imacs_i = SpectralElement(Empirical1D, points=imacs_i[:,0], lookup_table=imacs_i[:,1])
bp_imacs_z = SpectralElement(Empirical1D, points=imacs_z[:,0], lookup_table=imacs_z[:,1])

#obs_ig = Observation(LBG_with_IGM, bp_imacs_g, force = 'extrap')#, binset='array-like')
obs_ir = Observation(LBG, bp_imacs_r, force = 'extrap')#, binset='array-like')
obs_ii = Observation(LBG, bp_imacs_i, force = 'extrap')#,  binset='array-like')
obs_iz = Observation(LBG, bp_imacs_z, force = 'extrap')#,  binset='array-like')


#z=3.627, 3.912

print('Subaru HSC filters mag')
print('redshift', z)
print('r-band = ',obs_hr.effstim('abmag'))
print('i-band = ',obs_hi.effstim('abmag'))
print('z-band = ',obs_hz.effstim('abmag'))
print('IMACS filters mag')
print('redshift', z)
print('r-band = ',obs_ir.effstim('abmag'))
print('i-band = ',obs_ii.effstim('abmag'))
print('z-band = ',obs_iz.effstim('abmag'))

print('color delve')
print('ri-color = ',obs_ir.effstim('abmag')-obs_ii.effstim('abmag'))
print('color  Subaru HSC')
print('ri-color = ',obs_hr.effstim('abmag')-obs_hi.effstim('abmag'))

print('color delve ')
print('iz-color = ',obs_ii.effstim('abmag')-obs_iz.effstim('abmag'))

print('color  Subaru HSC')
print('iz-color = ',obs_hi.effstim('abmag')-obs_hz.effstim('abmag'))