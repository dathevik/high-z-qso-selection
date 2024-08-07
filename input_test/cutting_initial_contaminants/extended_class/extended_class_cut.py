import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Read data
data_qso_araa = ascii.read('delve_class_z_araa.dat')
data_qso_milliquas = ascii.read('delve_class_z_milliquas.dat')
data_qso_yang = ascii.read('delve_class_z_yang.dat')
data_bd = ascii.read('delve_class_z_bd.dat')

ext_z_yang = data_qso_yang["extended_class_z"]
ext_z_milliquas = data_qso_milliquas["extended_class_z"]
ext_z_araa = data_qso_araa["extended_class_z"]
ext_z_bd = data_bd["extended_class_z"]

# Custom bins
bins = [-10, -9, -8, 0, 1, 2, 3, 4]  # Include 4 to capture the last bin edge at 3
plt.figure(figsize=(10,6), dpi=200 )
# Create histogram
plt.hist(ext_z_bd, bins=bins, edgecolor='gray', align='left', color='white', alpha=0.3,histtype='stepfilled',  stacked='True', hatch='||', label='Brown Dwarfs')
plt.hist(ext_z_milliquas, bins=bins, align='left', color='white', alpha=0.3, histtype='stepfilled', edgecolor='green', stacked='True',hatch='//', label='Quasars Flesch+2023')
plt.hist(ext_z_yang, bins=bins, edgecolor='blue', align='left', color='white', alpha=0.3, histtype='stepfilled',  stacked='True',hatch='\\\\', label='Quasars Yang+2023')
plt.hist(ext_z_araa, bins=bins, edgecolor='red', align='left', color='white', alpha=0.3,histtype='stepfilled',  stacked='True', hatch='----', label='Quasars Fan+2022')
# Set x-ticks to match the bin centers
plt.xticks(bins)
plt.legend()
# Set x-ticks to match the bin centers
plt.xticks([-9, 0, 1, 2, 3])
plt.xlim(-10, 4)

# Adding titles and labels
plt.xlabel('Extended class z')
plt.ylabel('Frequency')
plt.legend(loc="upper left")

# Show the plot
# plt.show()
plt.savefig("QSO_BD_extended_z.png")