import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import matplotlib.gridspec as gridspec

# Read data
data_qso_araa = ascii.read('araa_decals.dat')
data_qso_milliquas = ascii.read('milliquas_decals.dat')
data_qso_yang = ascii.read('yang_decals.dat')
data_bd = ascii.read('bd_decals.dat', encoding='latin-1')

g_mag_yang = data_qso_yang["mag_g"]
g_mag_milliquas = data_qso_milliquas["mag_g"]
g_mag_araa = data_qso_araa["mag_g"]
g_mag_bd = data_bd["mag_g"]

g_snr_yang = data_qso_yang["snr_g"]
g_snr_milliquas = data_qso_milliquas["snr_g"]
g_snr_araa = data_qso_araa["snr_g"]
g_snr_bd = data_bd["snr_g"]

r_mag_yang = data_qso_yang["mag_r"]
r_mag_milliquas = data_qso_milliquas["mag_r"]
r_mag_araa = data_qso_araa["mag_r"]
r_mag_bd = data_bd["mag_r"]

i_mag_yang = data_qso_yang["mag_i"]
i_mag_milliquas = data_qso_milliquas["mag_i"]
i_mag_araa = data_qso_araa["mag_i"]
i_mag_bd = data_bd["mag_i"]

z_mag_yang = data_qso_yang["mag_z"]
z_mag_milliquas = data_qso_milliquas["mag_z"]
z_mag_araa = data_qso_araa["mag_z"]
z_mag_bd = data_bd["mag_z"]

redshift_yang = data_qso_yang["redshift"]
redshift_milliquas = data_qso_milliquas["Z"]
redshift_araa = data_qso_araa["redshift"]

# Read DELVE data for mag_g
data_qso_araa_delve = ascii.read('../g_band_detection/delve_araa_g_band.dat')
data_qso_milliquas_delve = ascii.read('../g_band_detection/delve_miliquas_g_band.dat')
data_qso_yang_delve = ascii.read('../g_band_detection/delve_yang_g_band.dat')
data_bd_delve = ascii.read('../g_band_detection/delve_bd_g_band.dat')

g_mag_yang_delve = data_qso_yang_delve["mag_auto_g"]
g_mag_milliquas_delve = data_qso_milliquas_delve["mag_auto_g"]
g_mag_araa_delve = data_qso_araa_delve["mag_auto_g"]
g_mag_bd_delve = data_bd_delve["mag_auto_g"]

# Function to find matching indices based on RA and Dec
def find_matches(ra1, dec1, ra2, dec2, tolerance=0.0003):  # tolerance ~1 arcsec
    matches1 = []
    matches2 = []
    for i, (r1, d1) in enumerate(zip(ra1, dec1)):
        # Find all points within tolerance
        matches = np.where((np.abs(ra2 - r1) < tolerance) & (np.abs(dec2 - d1) < tolerance))[0]
        if len(matches) > 0:
            # Take the closest match if multiple matches found
            dists = np.sqrt((ra2[matches] - r1)**2 + (dec2[matches] - d1)**2)
            closest = matches[np.argmin(dists)]
            matches1.append(i)
            matches2.append(closest)
    return np.array(matches1), np.array(matches2)

# Find matching indices for each catalog
bd_match, bd_delve_match = find_matches(data_bd['ra'], data_bd['dec'], 
                                      data_bd_delve['ra_delve'], data_bd_delve['dec_delve'])
milliquas_match, milliquas_delve_match = find_matches(data_qso_milliquas['ra'], data_qso_milliquas['dec'],
                                                     data_qso_milliquas_delve['ra_delve'], data_qso_milliquas_delve['dec_delve'])
yang_match, yang_delve_match = find_matches(data_qso_yang['ra'], data_qso_yang['dec'],
                                          data_qso_yang_delve['ra_delve'], data_qso_yang_delve['dec_delve'])
araa_match, araa_delve_match = find_matches(data_qso_araa['ra'], data_qso_araa['dec'],
                                          data_qso_araa_delve['ra_delve'], data_qso_araa_delve['dec_delve'])

# Get matched magnitudes
g_mag_bd = data_bd['mag_g'][bd_match]
g_mag_bd_delve = data_bd_delve['mag_auto_g'][bd_delve_match]
g_mag_milliquas = data_qso_milliquas['mag_g'][milliquas_match]
g_mag_milliquas_delve = data_qso_milliquas_delve['mag_auto_g'][milliquas_delve_match]
g_mag_yang = data_qso_yang['mag_g'][yang_match]
g_mag_yang_delve = data_qso_yang_delve['mag_auto_g'][yang_delve_match]
g_mag_araa = data_qso_araa['mag_g'][araa_match]
g_mag_araa_delve = data_qso_araa_delve['mag_auto_g'][araa_delve_match]

# Print number of matches found
print(f"Number of matches found:")
print(f"Brown Dwarfs: {len(bd_match)}")
print(f"Milliquas: {len(milliquas_match)}")
print(f"Yang: {len(yang_match)}")
print(f"ARAA: {len(araa_match)}")

# Filter out non-matching entries (where either catalog has mag=99 or mag=0)
valid_bd = (g_mag_bd != 99) & (g_mag_bd != 0) & (g_mag_bd_delve != 99) & (g_mag_bd_delve != 0)
valid_milliquas = (g_mag_milliquas != 99) & (g_mag_milliquas != 0) & (g_mag_milliquas_delve != 99) & (g_mag_milliquas_delve != 0)
valid_yang = (g_mag_yang != 99) & (g_mag_yang != 0) & (g_mag_yang_delve != 99) & (g_mag_yang_delve != 0)
valid_araa = (g_mag_araa != 99) & (g_mag_araa != 0) & (g_mag_araa_delve != 99) & (g_mag_araa_delve != 0)

# Print number of valid matches
print(f"\nNumber of valid matches (excluding mag=99 or mag=0):")
print(f"Brown Dwarfs: {np.sum(valid_bd)}")
print(f"Milliquas: {np.sum(valid_milliquas)}")
print(f"Yang: {np.sum(valid_yang)}")
print(f"ARAA: {np.sum(valid_araa)}")

# DECALS g_mag vs DELVE g_mag comparison
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(g_mag_bd[valid_bd], g_mag_bd_delve[valid_bd], label="Brown Dwarfs Sample", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(g_mag_milliquas[valid_milliquas], g_mag_milliquas_delve[valid_milliquas], label="Fl23", s=16, alpha=0.4, marker="x", color="blue")
plt.scatter(g_mag_yang[valid_yang], g_mag_yang_delve[valid_yang], label="Y23", s=16, alpha=0.4, marker="x", color="red")
plt.scatter(g_mag_araa[valid_araa], g_mag_araa_delve[valid_araa], label="F23", s=16, alpha=0.4, marker="x", color="green")
plt.xlabel("DECALS g_mag")
plt.ylabel("DELVE g_mag")
plt.legend(loc="upper left")
plt.xlim(20, 32)
plt.ylim(20, 32)
# Add diagonal line
plt.plot([20, 32], [20, 32], 'k--', alpha=0.5)
plt.show()
#plt.savefig("BD_QSO_mag_g_decals_delve.png")
plt.close()

# Create plot with g_mag vs redshift
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_qso_milliquas['Z'], data_qso_milliquas['mag_g'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['redshift'], data_qso_yang['mag_g'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['redshift'], data_qso_araa['mag_g'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.xlabel("Redshift")
plt.ylabel("DECALS g_mag")
plt.legend(loc="upper left")
plt.ylim(20, 32)
plt.xlim(4, 7)
plt.grid(True, alpha=0.3)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tight_layout()
#plt.show()
plt.close()

# G-I vs redshift
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_qso_milliquas['Z'], data_qso_milliquas['mag_g'] - data_qso_milliquas['mag_i'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['redshift'], data_qso_yang['mag_g'] - data_qso_yang['mag_i'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['redshift'], data_qso_araa['mag_g'] - data_qso_araa['mag_i'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.xlabel("Redshift")
plt.ylabel("g-i")
plt.legend(loc="upper left")
plt.xlim(4, 7)
plt.ylim(-2, 12)
# plt.show()
#plt.savefig("BD_QSO_g_i_decals.png")
plt.close()

# G-R vs redshift
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_qso_milliquas['Z'], data_qso_milliquas['mag_g'] - data_qso_milliquas['mag_r'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['redshift'], data_qso_yang['mag_g'] - data_qso_yang['mag_r'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['redshift'], data_qso_araa['mag_g'] - data_qso_araa['mag_r'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.xlabel("Redshift")
plt.ylabel("g-r")
plt.legend(loc="upper left")
plt.xlim(4, 7)
plt.ylim(-2, 12)
# plt.show()
#plt.savefig("BD_QSO_g_r_decals.png")
plt.close()

# G-Z vs redshift
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_qso_milliquas['Z'], data_qso_milliquas['mag_g'] - data_qso_milliquas['mag_z'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['redshift'], data_qso_yang['mag_g'] - data_qso_yang['mag_z'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['redshift'], data_qso_araa['mag_g'] - data_qso_araa['mag_z'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.xlabel("Redshift")
plt.ylabel("g-z")
plt.legend(loc="upper left")
plt.xlim(4, 7)
plt.ylim(-2, 12)
# plt.show()
#plt.savefig("BD_QSO_g_z_decals.png")
plt.close()

# G-I vs G-R
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_bd['mag_g'] - data_bd['mag_i'], data_bd['mag_g'] - data_bd['mag_z'], label="Brown Dwarfs Sample", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(data_qso_milliquas['mag_g'] - data_qso_milliquas['mag_i'], data_qso_milliquas['mag_g'] - data_qso_milliquas['mag_z'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['mag_g'] - data_qso_yang['mag_i'], data_qso_yang['mag_g'] - data_qso_yang['mag_z'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['mag_g'] - data_qso_araa['mag_i'], data_qso_araa['mag_g'] - data_qso_araa['mag_z'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.xlabel("g-i")
plt.ylabel("g-z")
plt.legend(loc="upper left")
plt.xlim(0, 12)
plt.ylim(0, 12)
# plt.show()
#plt.savefig("BD_QSO_g_z_g_i_decals.png")
plt.close()

# G-I vs snr_g
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_bd['mag_g'] - data_bd['mag_i'], data_bd['snr_g'], label="Brown Dwarfs Sample", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(data_qso_milliquas['mag_g'] - data_qso_milliquas['mag_i'], data_qso_milliquas['snr_g'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['mag_g'] - data_qso_yang['mag_i'], data_qso_yang['snr_g'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['mag_g'] - data_qso_araa['mag_i'], data_qso_araa['snr_g'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.xlabel("g-i")
plt.ylabel("SNR_g")
plt.axhline(y=3, color='gray', linestyle='--', label='SNR < 3 cut')
plt.axvline(x=2, color='red', linestyle='--', label='SNR >3 and g-i > 2 cut')
plt.legend(loc="upper left")
plt.xlim(0, 12)
plt.ylim(0, 25)
# plt.show()
#plt.savefig("BD_QSO_snr_g_i_cut.png")
plt.close()

# G-Z vs snr_g
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_bd['mag_g'] - data_bd['mag_z'], data_bd['snr_g'], label="Brown Dwarfs Sample", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(data_qso_milliquas['mag_g'] - data_qso_milliquas['mag_z'], data_qso_milliquas['snr_g'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['mag_g'] - data_qso_yang['mag_z'], data_qso_yang['snr_g'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['mag_g'] - data_qso_araa['mag_z'], data_qso_araa['snr_g'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.axhline(y=3, color='gray', linestyle='--', label='SNR < 3 cut')
plt.axvline(x=2, color='red', linestyle='--', label='SNR >3 and g-z > 2 cut')
plt.xlabel("g-z")
plt.ylabel("SNR_g")
plt.legend(loc="upper left")
plt.xlim(0, 12)
plt.ylim(0, 25)
# plt.show()
#plt.savefig("BD_QSO_snr_g_z_cut.png")
plt.close()

# DECALS g_mag vs SNR_g
plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(data_bd['mag_g'], data_bd['snr_g'], label="Brown Dwarfs Sample", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(data_qso_milliquas['mag_g'], data_qso_milliquas['snr_g'], label="Fl23", s=8, alpha=0.4, marker="o", color="blue")
plt.scatter(data_qso_yang['mag_g'], data_qso_yang['snr_g'], label="Y23", s=8, alpha=0.4, marker="o", color="red")
plt.scatter(data_qso_araa['mag_g'], data_qso_araa['snr_g'], label="F23", s=8, alpha=0.4, marker="o", color="green")
plt.axhline(y=3, color='gray', linestyle='--', label='SNR < 3 cut')
plt.xlabel("DECaLS g_mag")
plt.ylabel("SNR_g")
plt.legend(loc="upper right")
plt.xlim(20, 32)
plt.ylim(0, 50)
plt.grid(True, alpha=0.3)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tight_layout()
#plt.show()
#plt.savefig("BD_QSO_mag_g_snr_g_decals.png")
plt.close()