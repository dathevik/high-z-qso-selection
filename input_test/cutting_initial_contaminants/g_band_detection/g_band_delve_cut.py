import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Read data
data_qso_araa = ascii.read('delve_araa_g_band.dat')
data_qso_milliquas = ascii.read('delve_miliquas_g_band.dat')
data_qso_yang = ascii.read('delve_yang_g_band.dat')
data_bd = ascii.read('delve_bd_g_band.dat')

g_mag_yang = data_qso_yang["mag_auto_g"]
g_mag_milliquas = data_qso_milliquas["mag_auto_g"]
g_mag_araa = data_qso_araa["mag_auto_g"]
g_mag_bd = data_bd["mag_auto_g"]

g_magerr_yang = data_qso_yang["magerr_auto_g"]
g_magerr_milliquas = data_qso_milliquas["magerr_auto_g"]
g_magerr_araa = data_qso_araa["magerr_auto_g"]
g_magerr_bd = data_bd["magerr_auto_g"]

snr_yang = 1.0857/g_magerr_yang
snr_milliquas = 1.0857/g_magerr_milliquas
snr_araa = 1.0857/g_magerr_araa
snr_bd = 1.0857/g_magerr_bd

plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(g_mag_bd,snr_bd, label="Brown Dwarfs", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(g_mag_milliquas, snr_milliquas, label="Quasars Flesch+2023", s=16, alpha=0.4, marker="x", color="blue")
plt.scatter(g_mag_araa, snr_araa, label="Quasars Fan+2022", s=50, alpha=0.6, marker="x", color="green")
plt.scatter(g_mag_yang, snr_yang, label="Quasars Yang+2023", s=50, alpha=0.6, marker="x", color="red")
plt.xlabel("DELVE g_mag")
plt.ylabel("Signal to Noise")
plt.legend(loc="upper left")
plt.xlim(18, 30)
plt.ylim(0, 100)
# plt.show()
plt.savefig("BD_QSO_mag_g_snr.png")
plt.close()

bins = [16, 18, 20, 22, 24, 26, 30, 97, 99, 100]  # Include 4 to capture the last bin edge at 3
plt.figure(figsize=(10,6), dpi=200 )
# Create histogram
plt.hist(g_mag_bd, bins=bins, edgecolor='gray', align='left', color='white', alpha=0.5,histtype='stepfilled',  stacked='True', hatch='||', label='Brown Dwarfs')
plt.hist(g_mag_milliquas, bins=bins, align='left', color='white', alpha=0.3, histtype='stepfilled', edgecolor='green', stacked='True',hatch='//', label='Quasars Flesch+2023')
plt.hist(g_mag_yang, bins=bins, edgecolor='blue', align='left', color='white', alpha=0.3, histtype='stepfilled',  stacked='True',hatch='\\\\', label='Quasars Yang+2023')
plt.hist(g_mag_araa, bins=bins, edgecolor='red', align='left', color='white', alpha=0.3,histtype='stepfilled',  stacked='True', hatch='----', label='Quasars Fan+2022')
# Set x-ticks to match the bin centers
plt.xticks([15, 20, 25, 99])
plt.xlabel("DELVE g_mag")
plt.legend()
plt.xlim(18, 100)
plt.savefig("BD_QSO_mag_g_hist.png")
# plt.show()