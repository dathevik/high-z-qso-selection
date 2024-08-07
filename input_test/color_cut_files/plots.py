import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.units as u

#Convert fluxes to magnitudes
def flux_to_mag(flux, flux_zero):
    return (-2.5 * np.log10(flux / flux_zero)) * u.ABmag

# Zero-point flux values in mJy
flux_zero_r = 3.631
flux_zero_I = 3.631
flux_zero_z = 3.631

#Reading files and preparing plotting parameters in a function
def plot_template(file_name):
    template_path = os.path.abspath(f'{file_name}')
    data_temp = ascii.read(template_path)
    print('------------------ TEMPLATE FILE -------------------\n', data_temp)
    r_decam_flux = data_temp.columns[2]
    i_decam_flux = data_temp.columns[3]
    z_decam_flux = data_temp.columns[4]
    redshift = data_temp.columns[0]
    z_decam_mag = []
    r_decam_mag = []
    i_decam_mag = []

    x_mag = []
    y_mag = []

    for mag in range(len(z_decam_flux)):
        mag_z = z_decam_flux[mag]
        mag_r = r_decam_flux[mag]
        mag_i = i_decam_flux[mag]
        z_decam_mag.append(mag_z)
        r_decam_mag.append(mag_r)
        i_decam_mag.append(mag_i)
        r_i = mag_r - mag_i
        i_z = mag_i - mag_z
        x_mag.append(r_i)
        y_mag.append(i_z)
    print('------------------------------------------------------------')
    print('---magnitudes of r band--- \n', r_decam_mag)
    print('------------------------------------------------------------')
    print('---magnitudes of i band--- \n', i_decam_mag)
    print('------------------------------------------------------------')
    print('---magnitudes of z band--- \n', z_decam_mag)
    print('------------------------------------------------------------')
    print('---r-i color difference--- \n', x_mag)
    print('------------------------------------------------------------')
    print('---i-z color difference--- \n', y_mag)
    print('------------------------------------------------------------')
    print('Max magnitude z', max(z_decam_mag))
    print('Max magnitude r', max(r_decam_mag))
    print('Max magnitude i', max(i_decam_mag))
# Plotting the template points and corresponding redshifts as a subplot
    plt.plot(x_mag, y_mag, label=f'{file_name}')
    plt.scatter(x_mag, y_mag, color='black', s=4)
    for i, txt in enumerate(redshift):
        if y_mag[i] > 1.3:
            plt.annotate(txt, (x_mag[i], y_mag[i]), fontsize=4)
    #plt.savefig(f'{save_name}template_mag_color_difference.png')
# fig, ax = plt.subplots()
# rectangle = plt.Rectangle((1, 0), 2, 2, fill='False', edgecolor='black', facecolor='white', linestyle='dashed')
# ax.add_patch(rectangle)
#plotting quasars datapoints
araa_path = os.path.abspath('araa_delve_query.dat')
data_araa = ascii.read(araa_path)
r_decam_mag_qso = data_araa.columns[0]
i_decam_mag_qso = data_araa.columns[1]
z_decam_mag_qso = data_araa.columns[2]
x_mag_qso = r_decam_mag_qso - i_decam_mag_qso
y_mag_qso = i_decam_mag_qso - z_decam_mag_qso
redshift_qso = data_araa.columns[9]
# for i, txt in enumerate(redshift_qso):
#      ax.annotate(txt, (x_mag_qso[i], y_mag_qso[i]), fontsize=3.5)

#plotting brown dwarf datapoints
bd_path = os.path.abspath('bd_delve_query.dat')
data_bd = ascii.read(bd_path)
print('------------------ OBJECT FILE -------------------\n', data_bd)
r_decam_mag_bd = data_bd.columns[6]
i_decam_mag_bd = data_bd.columns[7]
z_decam_mag_bd = data_bd.columns[8]
x_mag_bd = r_decam_mag_bd - i_decam_mag_bd
y_mag_bd = i_decam_mag_bd - z_decam_mag_bd
print('------------------ Brown Dwarfs list r-i and i-z -------------------\n')
print('r-i', x_mag_bd)
print('i-z is', y_mag_bd)

# plotting templates (Selsing)
sels_path = os.path.abspath('Selsing2015_fluxes_mJy_202301.dat')
data_sels = ascii.read(sels_path)
redshift_sels = data_sels.columns[1]
r_decam_flux_sels = data_sels.columns[3]
i_decam_flux_sels = data_sels.columns[4]
z_decam_flux_sels = data_sels.columns[5]
r_decam_mag_sels = flux_to_mag(r_decam_flux_sels, flux_zero_r)
i_decam_mag_sels = flux_to_mag(i_decam_flux_sels, flux_zero_I)
z_decam_mag_sels = flux_to_mag(z_decam_flux_sels, flux_zero_z)
x_mag_sels = r_decam_mag_sels - i_decam_mag_sels
y_mag_sels = i_decam_mag_sels - z_decam_mag_sels


# plotting BD templates (Roberto's sample)
roberto_path = os.path.abspath('BDRA_fluxes_mJy_202301.dat')
data_roberto = ascii.read(roberto_path)
# redshift_sels = data_sels.columns[1]
r_decam_flux_rob = data_roberto.columns[2]
i_decam_flux_rob = data_roberto.columns[3]
z_decam_flux_rob = data_roberto.columns[4]
r_decam_mag_rob = flux_to_mag(r_decam_flux_rob, flux_zero_r)
i_decam_mag_rob = flux_to_mag(i_decam_flux_rob, flux_zero_I)
z_decam_mag_rob = flux_to_mag(z_decam_flux_rob, flux_zero_z)
x_mag_rob = r_decam_mag_rob - i_decam_mag_rob
y_mag_rob = i_decam_mag_rob - z_decam_mag_rob
print('------------------ Brown Dwarf template r-i and i-z -------------------\n')
print('r-i', x_mag_rob)
print('i-z is', y_mag_rob)

# plotting everything together as subplots
plot_template('qsogen_mags_emlines2.0_ebv+025.dat')
plot_template('qsogen_mags_emlines2.0_ebv+000.dat')
plot_template('qsogen_mags_emlines-2.0_ebv+025.dat')
plot_template('qsogen_mags_emlines-2.0_ebv+000.dat')
plt.xlabel('r-i')
plt.ylabel('i-z')
# ax.plot(x_mag_sels, y_mag_sels, color='darkgreen', label='Selsing templates')
# ax.plot(x_mag_rob, y_mag_rob, color='darkblue', label='BD templates')
plt.xlim(-1, 3)
# ax.scatter(x_mag_bd, y_mag_bd, s=8, color='saddlebrown', label='known brown dwarfs')
# ax.scatter(x_mag_qso, y_mag_qso, s=8, color='khaki', label='known quasars')
plt.legend(fontsize='7')
plt.savefig('plots/all_temp+bd+qso_datapoints_in_one.png', dpi=400)

