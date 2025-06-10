import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Read data
data_qso_araa = ascii.read('delve_class_z_araa.dat')
data_qso_milliquas = ascii.read('delve_class_z_milliquas.dat')
data_qso_yang = ascii.read('delve_class_z_yang.dat')
data_bd = ascii.read('delve_class_z_bd.dat')
data_sdss = ascii.read('delve_class_z_sdss12.dat')

ext_z_yang = data_qso_yang["extended_class_z"]
ext_z_milliquas = data_qso_milliquas["extended_class_z"]
ext_z_araa = data_qso_araa["extended_class_z"]
ext_z_bd = data_bd["extended_class_z"]
ext_z_sdss = data_sdss["extended_class_z"]

# Custom bins
bins = [-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
plt.figure(figsize=(12, 6), dpi=300)
plt.grid(True, linestyle='--', alpha=0.3)

# Create histogram
plt.hist(ext_z_bd, bins=bins, edgecolor='gray', align='left', color='grey', alpha=0.3, histtype='stepfilled', stacked='True', hatch='||', label='BD Sample')
plt.hist(ext_z_milliquas, bins=bins, align='left', color='lightgreen', alpha=0.7, histtype='stepfilled', edgecolor='darkgreen', stacked='True', hatch='//', label='Fl23 Sample')
plt.hist(ext_z_yang, bins=bins, edgecolor='blue', align='left', color='lightblue', alpha=0.7, histtype='stepfilled', stacked='True', hatch='\\\\', label='Y23 Sample')
plt.hist(ext_z_araa, bins=bins, edgecolor='red', align='left', color='pink', alpha=0.4, histtype='stepfilled', stacked='True', hatch='----', label='F23 Sample')
plt.hist(ext_z_sdss, bins=bins, edgecolor='purple', align='left', color='lavender', alpha=0.7, histtype='stepfilled', stacked='True', hatch='++', label='Galaxies')

plt.xticks([-1, 0, 1, 2, 3], fontsize=18)
plt.xlim(-1, 4)
plt.ylim(0, 4000)

plt.yticks([0, 1000, 2000, 3000, 4000], fontsize=18)

plt.xlabel('extended_class_z_delve', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.legend(loc="upper right", fontsize=18, framealpha=0.8)

plt.tight_layout()
plt.savefig("QSO_BD_extended_z.png")