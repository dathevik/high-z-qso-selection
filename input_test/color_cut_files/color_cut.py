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
def plot_template(file_name, labels):
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
    #plt.plot(x_mag, y_mag, label=labels)
    # plt.scatter(x_mag, y_mag, color='black', s=4)
    #     for i, txt in enumerate(redshift):
    #         if y_mag[i] > 1:
    #             plt.annotate(txt, (x_mag[i], y_mag[i]), fontsize=8)


def plot_template_redshift(file_list, type, labels):
    ri_min = []
    ri_max = []
    iz_min = []
    iz_max = []
    if type == 0:
        for file, label in zip(file_list, labels):
            template_path = os.path.abspath(f'{file}')
            data_temp = ascii.read(template_path)
            print('------------------ TEMPLATE FILE -------------------\n', data_temp)
            r_decam_flux = data_temp.columns[2]
            i_decam_flux = data_temp.columns[3]
            redshift = data_temp.columns[0]
            ri = r_decam_flux - i_decam_flux
            ri_mags = []
            for mag in range(len(redshift[:120])):
                if ri[mag] > 1.3:
                    ri_mags.append(ri[mag])
                    ri_min.append(min(ri_mags))
                    ri_max.append(max(ri_mags))
            plt.plot(redshift, ri, label=label)
            plt.legend(fontsize='5')
        plt.axvline(x=6, color='black', linestyle='dashed',  label='Redshift range')
        print('the minimum ri of all templates is', min(ri_min))
        print('the maximum ri of all templates is', max(ri_max))
        plt.axhline(y=min(ri_min), color='red', linestyle='dashed', label='4MOST threshold')
        # plt.axhline(y=min(ri_min)-0.152, color='blue', linestyle='dashed', label='COSMOS threshold')
        plt.axvline(x=4.9, color='black', linestyle='dashed')
        plt.xlabel('redshift', fontsize=14)
        plt.ylabel('r-i', fontsize=14)
        plt.xticks([0, 1,  2,  3, 4, 4.9, 6, 7, 8], fontsize='9')
        plt.yticks([-1, 0, 1, min(ri_min), 2, 3, 4], fontsize='9') # for cosmos  min(ri_min)-0.152
        plt.ylim(-1, 4)
        plt.xlim(0, 8)

    elif type == 1:
        for file, label in zip(file_list, labels):
            template_path = os.path.abspath(f'{file}')
            data_temp = ascii.read(template_path)
            print('------------------ TEMPLATE FILE -------------------\n', data_temp)
            i_decam_flux = data_temp.columns[3]
            z_decam_flux = data_temp.columns[4]
            redshift = data_temp.columns[0]
            iz = i_decam_flux - z_decam_flux
            iz_mags = []
            for mag in range(len(redshift[:130])):
                if iz[mag] > 1.5:
                    iz_mags.append(iz[mag])
                    iz_min.append(min(iz_mags))
                    iz_max.append(max(iz_mags))
            plt.plot(redshift, iz, label=label)
            plt.legend(fontsize='5')
        print('the minimum iz of all templates is', min(iz_min))
        print('the maximum iz of all templates is', max(iz_max))
        plt.axhline(y=min(iz_min), color='red', linestyle='dashed', label='4MOST threshold')
        # plt.axhline(y=min(iz_min)+0.106, color='blue', linestyle='dashed', label='COSMOS threshold')
        plt.axvline(x=6, color='black', linestyle='dashed', label='Redshift range')
        plt.axvline(x=6.97, color='black', linestyle='dashed')
        plt.xlabel('redshift', fontsize=14)
        plt.ylabel('i-z', fontsize=14)
        plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8], fontsize='8')
        plt.yticks([-1, 0, 1, min(iz_min), 2, 3, 4], fontsize='8') # for cosmos min(iz_min)+0.106
        plt.ylim(-1, 4)
        plt.xlim(0, 8)

#plotting quasars datapoints
araa_path = os.path.abspath('araa_delve_query.dat')
data_araa = ascii.read(araa_path)
r_mag_qso = data_araa.columns[0]
i_mag_qso = data_araa.columns[1]
z_mag_qso = data_araa.columns[2]
ri_mag_qso = r_mag_qso - i_mag_qso
iz_mag_qso = i_mag_qso - z_mag_qso
redshift_qso = data_araa.columns[9]

# #plotting cosmos2020 datapoints
# cosmos_path = os.path.abspath('cosmos_stars_2020_griz.dat')
# data_cosmos = ascii.read(cosmos_path)
# r_subaru_cosmos = data_cosmos.columns[4]
# i_subaru_cosmos = data_cosmos.columns[5]
# z_subaru_cosmos = data_cosmos.columns[6]
# ri_mag_cosmos = r_subaru_cosmos - i_subaru_cosmos
# iz_mag_cosmos = i_subaru_cosmos - z_subaru_cosmos
# redshift_cosmos = data_cosmos.columns[7]

#plotting brown dwarf datapoints
bd_path = os.path.abspath('bd_delve_query.dat')
data_bd = ascii.read(bd_path)
print('------------------ OBJECT FILE -------------------\n', data_bd)
r_decam_mag_bd = data_bd.columns[6]
i_decam_mag_bd = data_bd.columns[7]
z_decam_mag_bd = data_bd.columns[8]
ri_mag_bd = r_decam_mag_bd - i_decam_mag_bd
iz_mag_bd = i_decam_mag_bd - z_decam_mag_bd
print('------------------ Brown Dwarfs list r-i and i-z -------------------\n')
print('r-i', ri_mag_bd)
print('i-z is',iz_mag_bd)

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
ri_sels = r_decam_mag_sels - i_decam_mag_sels
iz_sels = i_decam_mag_sels - z_decam_mag_sels


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
ri_rob = r_decam_mag_rob - i_decam_mag_rob
iz_rob = i_decam_mag_rob - z_decam_mag_rob
print('------------------ Brown Dwarf template r-i and i-z -------------------\n')
print('r-i is', ri_rob)
print('i-z is', iz_rob)

# plotting everything together
# plot_template('qsogen_mags_emlines2.0_ebv+025.dat')
# plot_template('qsogen_mags_emlines2.0_ebv+000.dat')
# plot_template('qsogen_mags_emlines-2.0_ebv+025.dat')
# plot_template('qsogen_mags_emlines-2.0_ebv+000.dat')
# plt.xlabel('r-i')
# plt.ylabel('i-z')
# # plt.plot(ri_sels, iz_sels, color='darkgreen', label='Selsing templates')
# # plt.plot(ri_rob, iz_rob, color='darkblue', label='BD templates')
# plt.xlim(-1, 5)
# plt.ylim(-1, 6)
# plt.scatter(ri_mag_cosmos, iz_mag_cosmos, s=0.01, color='black', label='cosmos sources')
# # plt.scatter(ri_mag_bd, iz_mag_bd, s=10, color='lime', label='known brown dwarfs')
# # plt.scatter(ri_mag_qso, iz_mag_qso, s=10, color='yellow', label='known quasars')
# plt.legend(fontsize='7')
# plt.savefig('plots/all_temp+bd+qso_datapoints_cosmos.png', dpi=400)
# plt.close()


# plot ri or iz with redshift
# plt.scatter(redshift_cosmos, iz_mag_cosmos, s=0.01, color='gray', label='Cosmos 2020 sources')
matt_filelist = ['qsogen_mags_emlines2.0_ebv+025.dat', 'qsogen_mags_emlines-2.0_ebv+010.dat', 'qsogen_mags_emlines0.0_ebv+000.dat']
label_temp = ['qsogen template, emlines: 2, ebv: 0.25', 'qsogen template, emlines: -2, ebv: 0.1', 'qsogen template, emlines: 0, ebv: 0']
plt.scatter(redshift_qso, ri_mag_qso, s=8, color='purple', label='known quasars')

# for i, txt in enumerate(redshift_qso):
#     if i % 3==0:
#         plt.annotate(txt, (redshift_qso[i], ri_mag_qso[i]), fontsize=8, color='black')
plt.plot(redshift_sels, ri_sels, label='Selsing template')
plot_template_redshift(matt_filelist, 0, label_temp)
plt.legend(fontsize='9')
plt.show()
# plt.savefig('plots/redshift_iz_numbers.png', dpi=400)
plt.close()

# print(iz_mag_cosmos, 'cosmos_iz')
# print(ri_mag_cosmos, 'cosmos_ri')
