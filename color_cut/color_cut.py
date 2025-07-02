import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.units as u

# Set the style for publication-quality plots
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['figure.figsize'] = [12, 10]
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

#Convert fluxes to magnitudes
def flux_to_mag(flux, flux_zero):
    return (-2.5 * np.log10(flux / flux_zero)) * u.ABmag

# Zero-point flux values in mJy
flux_zero_r = 3.631
flux_zero_I = 3.631
flux_zero_z = 3.631

# Define colorblind-friendly colors
template_colors = ['#0072B2', '#D55E00', '#009E73']  # Blue, Orange, Green for templates
selsing_color = '#CC79A7'  # Pink
qso_color = '#FFD700'      # Bright gold for Fan sample
fan_color = '#009E73'      # Green
cut_color = '#CC79A7'      # Pink
yang_color = '#4169E1'     # Royal Blue for Yang sample
milliquas_color = '#FF4500'  # Orange Red for Milliquas sample

# Create the main plot
fig, ax = plt.subplots(figsize=(12, 10))

# Plot QSO templates
matt_filelist = ['qsogen_mags_emlines2.0_ebv+025.dat', 'qsogen_mags_emlines-2.0_ebv+010.dat', 'qsogen_mags_emlines0.0_ebv+000.dat']
label_temp = ['QSO template, emlines: 2, E(B-V): 0.25 (Temple+2021)', 
              'QSO template, emlines: -2, E(B-V): 0.1 (Temple+2021)', 
              'QSO template, emlines: 0, E(B-V): 0 (Temple+2021)']

for file, label, color in zip(matt_filelist, label_temp, template_colors):
    template_path = os.path.abspath(f'{file}')
    data_temp = ascii.read(template_path)
    r_decam_flux = data_temp.columns[2]
    i_decam_flux = data_temp.columns[3]
    redshift = data_temp.columns[0]
    ri = r_decam_flux - i_decam_flux
    plt.plot(redshift, ri, label=label, color=color, linewidth=2, alpha=0.7)

# Plot Selsing templates
sels_path = os.path.abspath('Selsing2015_fluxes_mJy_202301.dat')
data_sels = ascii.read(sels_path)
redshift_sels = data_sels.columns[1]
r_decam_flux_sels = data_sels.columns[3]
i_decam_flux_sels = data_sels.columns[4]
r_decam_mag_sels = flux_to_mag(r_decam_flux_sels, flux_zero_r)
i_decam_mag_sels = flux_to_mag(i_decam_flux_sels, flux_zero_I)
ri_sels = r_decam_mag_sels - i_decam_mag_sels
plt.plot(redshift_sels, ri_sels, label='QSO template (Selsing+2016)', color=selsing_color, linewidth=2, alpha=0.7)

# Plot known quasars (Fan sample)
araa_path = os.path.abspath('araa_delve_query.dat')
data_araa = ascii.read(araa_path)
r_mag_qso = data_araa.columns[0]
i_mag_qso = data_araa.columns[1]
redshift_qso = data_araa.columns[9]
ri_mag_qso = r_mag_qso - i_mag_qso
plt.scatter(redshift_qso, ri_mag_qso, s=250, color=qso_color, alpha=0.7, label='F23 Sample', marker='*', edgecolor='black')

# Plot Yang 2023 quasar sample
yang_path = os.path.abspath('yang_results.dat')
try:
    data_yang = ascii.read(yang_path)
    r_mag_yang = data_yang['mag_auto_r']
    i_mag_yang = data_yang['mag_auto_i']
    redshift_yang = data_yang['QSO_z']
    ri_mag_yang = r_mag_yang - i_mag_yang
#    plt.scatter(redshift_yang, ri_mag_yang, s=60, color=yang_color, alpha=0.5, label='Y23 Sample', marker='s', edgecolor='black')
except Exception as e:
    print(f"Warning: Could not load Yang 2023 data. Error: {str(e)}")
    print("Available columns:", data_yang.colnames if 'data_yang' in locals() else "No data loaded")

# Plot Milliquas quasar sample
milliquas_path = os.path.abspath('milliquas_results.dat')
try:
    data_milliquas = ascii.read(milliquas_path)
    r_mag_milliquas = data_milliquas['mag_auto_r']
    i_mag_milliquas = data_milliquas['mag_auto_i']
    redshift_milliquas = data_milliquas['QSO_z']
    ri_mag_milliquas = r_mag_milliquas - i_mag_milliquas
#    plt.scatter(redshift_milliquas, ri_mag_milliquas, s=4, color=milliquas_color, alpha=0.4, label='Fl23 Sample', marker='o', edgecolor='black')
except Exception as e:
    print(f"Warning: Could not load Milliquas data. Error: {str(e)}")
    print("Available columns:", data_milliquas.colnames if 'data_milliquas' in locals() else "No data loaded")

# Plot Fan et al., 2023 quasar sample
# Assuming you have a file with Fan et al., 2023 data
fan_path = os.path.abspath('fan_et_al_2023_quasars.dat')
try:
    data_fan = ascii.read(fan_path)
    r_mag_fan = data_fan.columns[0]  # Adjust column indices as needed
    i_mag_fan = data_fan.columns[1]
    redshift_fan = data_fan.columns[2]
    ri_mag_fan = r_mag_fan - i_mag_fan
    plt.scatter(redshift_fan, ri_mag_fan, s=15, color=fan_color, alpha=0.6, label='F23 Sample', marker='o')
except:
    print("Warning: Could not load Fan et al., 2023 data. File may not exist or have different format.")

# Add color cut lines
# Replace the two separate lines with a rectangle for the selection region
plt.axvspan(4.5, 6.0, ymin=1.3/4, ymax=1.0, color=cut_color, alpha=0.2, label='Sources selected')

# Add dashed black frame around the selection region
plt.plot([4.5, 6.0, 6.0, 4.5, 4.5], [1.3, 1.3, 4.0, 4.0, 1.3], 'k--', linewidth=2)

# Add grid for better readability
ax.grid(True, linestyle='--', alpha=0.3)

# Improve axis labels
plt.xlabel('Redshift', fontsize=20)
plt.ylabel('r$_{delve}$ - i$_{delve}$', fontsize=20)

# Set plot limits
plt.xlim(4.5, 7)
plt.ylim(0, 4)

# Improve legend
plt.legend(fontsize=16, loc='upper right', frameon=True, framealpha=0.9, ncol=1)

# Add minor grid lines
ax.minorticks_on()
ax.grid(True, which='minor', linestyle=':', alpha=0.2)

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Show the plot without saving
#plt.show()
plt.savefig('plots/color_cut_fan.png', dpi=300)