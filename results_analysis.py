#----------- Making histogram from χ2 ratios calculated in sed_calculation.py 4MOST project script
#----------- Analysing results of the testing and improvements

from astropy.io import ascii, fits
import matplotlib.pyplot as plt
import os
import numpy as np

object_file = 'concat_F_and_BIC.fits'
output_path = os.path.abspath('/output_test')
objects_path = os.path.abspath('results_4most/F_and_BIC/' + object_file)
hdu = fits.open(objects_path)
data_results = hdu[1].data
ra = data_results.field(1)
dec = data_results.field(2)
BD_Chi2 = data_results.field(3)
QSO_Chi2 = data_results.field(5)
R_Chi2 = data_results.field(7)
F_test = data_results.field(8)
BIC = data_results.field(9)
print(data_results)

# ------- Scatter plot of ra and dec ratios
distr = plt.scatter(BIC, F_test,  edgecolor='black')


# ------- Histogram of F_test and BIC
# create new axes on the right and on the top of the current axes
dividerx = np.histogram(BIC)
dividery = np.histogram(F_test)
# below height and pad are in inches
ax_histx = dividerx.append_axes("top", 1.2, pad=0.1, sharex=BIC)
ax_histy = dividery.append_axes("right", 1.2, pad=0.1, sharey=F_test)
# make some labels invisible
ax_histx.xaxis.set_tick_params(labelbottom=False)
ax_histy.yaxis.set_tick_params(labelleft=False)

bins_x = np.arange(0.5, 3.5, 5.5)
bins_y = np.arange(1, 5, 10)
ax_histx.hist(np.log10(distr), bins=bins_x, color='grey', alpha=0.5)
ax_histy.hist(np.log10(distr), bins=bins_y, orientation='horizontal', alpha=0.5, color='grey')
plt.savefig(f"{output_path}/BIC_F_dist.png")