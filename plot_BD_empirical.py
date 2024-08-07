import os
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np


def flux_from_cgs_to_mJy(flux_cgs, ll):
    # ---- Obtain flux in mJy from flux in erg/s/cm2/A and considering effective wavelength of the filter in AA
    flux_mJy = (3.33564095E+04 * flux_cgs * ll ** 2) * 1E-13
    return flux_mJy


def clean_and_convert(value):
    # Remove any trailing characters that are not part of a number
    cleaned_value = value.strip().rstrip('.,')
    return float(cleaned_value)


def get_empirical_fluxes(empirical_fluxes):
    # Read the empirical flux file
    empirical_data = ascii.read(empirical_fluxes)
    spectral_types = empirical_data['Sp.T']

    # Convert the strings to floats, replacing commas with dots and cleaning the strings
    flux_columns = ['g', 'r', 'i', 'z', 'J', 'Y', 'H', 'K', 'W1', 'W2']
    fluxes = {}
    for col in flux_columns:
        if col in empirical_data.colnames:
            fluxes[col] = np.array([clean_and_convert(str(x).replace(',', '.')) for x in empirical_data[col]])
        else:
            fluxes[col] = None  # Column does not exist

    return spectral_types, fluxes


def plot_files_from_folder(folder_path, color_map, label, alpha, legend_color, empirical_fluxes):
    # Get the list of all files in the folder
    file_list = os.listdir(folder_path)

    # Filter the list to include only .txt files
    txt_files = [f for f in file_list if f.endswith('.txt')]

    # Create a color palette using a color map
    num_files = len(txt_files)
    colors = color_map(np.linspace(0, 1, num_files))

    # Get empirical fluxes
    spectral_types, fluxes = get_empirical_fluxes(empirical_fluxes)

    subtype_fluxes_and_colors = []

    # Plot all files with transparency and dashed lines using the color palette
    for i, file_name in enumerate(txt_files):
        file_path = os.path.join(folder_path, file_name)
        # Extract spectral type from filename
        spectral_type = file_name.split('.')[0]

        # Find the corresponding empirical fluxes
        idx = np.where(spectral_types == spectral_type)[0]
        if len(idx) > 0:
            idx = idx[0]

            # Read the file using astropy.io.ascii
            data = ascii.read(file_path)
            # Assuming the first and second columns are the ones to plot
            x_data = data.columns[0] * 10000
            y_data_cgs = data.columns[1]
            y_data = flux_from_cgs_to_mJy(y_data_cgs, x_data)

            # Plot the data with transparency, dashed lines, and color from the palette
            plt.plot(x_data, y_data, linestyle='--', alpha=alpha, color=colors[i])

            subtype_fluxes_and_colors.append((spectral_type, colors[i]))

    # Plot an invisible point to create a legend entry
    plt.plot([], [], linestyle='--', color=legend_color, alpha=alpha, label=label)

    return subtype_fluxes_and_colors


# Initialize a plot
plt.figure(figsize=(10, 6))

# Plotting all L-dwarfs with a specific alpha and legend color
folder_path_LDwarf = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/LDwarf/specs')
fluxes_LDwarf_empirical = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/L_type_wise_g_r.txt')
subtype_fluxes_and_colors_LDwarf = plot_files_from_folder(folder_path_LDwarf, plt.cm.Purples, 'LDwarf', alpha=0.3, legend_color='purple', empirical_fluxes=fluxes_LDwarf_empirical)

# Plotting all T-dwarfs with a different alpha and legend color
folder_path_TDwarf = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/TDwarf/specs')
fluxes_TDwarf_empirical = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/T_type_wise_g_r.txt')
subtype_fluxes_and_colors_TDwarf = plot_files_from_folder(folder_path_TDwarf, plt.cm.Oranges, 'TDwarf', alpha=0.3, legend_color='orange', empirical_fluxes=fluxes_TDwarf_empirical)

# Plotting all M-dwarfs with a different alpha and legend color
folder_path_MDwarf = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/MDwarf/specs')
fluxes_MDwarf_empirical = os.path.abspath('input_test/BD_templates_empirical_CM/BDspectra/M_type_wise_g_r.txt')
subtype_fluxes_and_colors_MDwarf = plot_files_from_folder(folder_path_MDwarf, plt.cm.Greens, 'MDwarf', alpha=0.3, legend_color='green', empirical_fluxes=fluxes_MDwarf_empirical)

# Function to plot empirical fluxes
def plot_empirical_fluxes(subtype_fluxes_and_colors, empirical_fluxes):
    # Read the empirical flux file
    empirical_data = ascii.read(empirical_fluxes)

    # Define central wavelengths
    central_wavelengths = {
        'g': 4798.3527009231575,  # g
        'r': 6407.493598028656,   # r
        'i': 7802.488114833454,   # i
        'z': 9144.625340022629,   # z
        'J': 12325.125694338809,  # J
        'Y': 10201.359507821942,  # Y
        'H': 16473.95843628733,   # H
        'Ks': 22045.772662096875,  # K
        'W1': 33791.878497259444,  # W1
        'W2': 46292.93969033106    # W2
    }

    # Loop through each subtype and plot its corresponding fluxes
    for subtype, color in subtype_fluxes_and_colors:
        # Extract fluxes for the current subtype from the empirical data
        subtype_fluxes = empirical_data[empirical_data['Sp.T'] == subtype]

        # If subtype data is found, plot the fluxes
        if len(subtype_fluxes) > 0:
            for flux_label, wavelength in central_wavelengths.items():
                if flux_label in subtype_fluxes.colnames:
                    flux_value = clean_and_convert(str(subtype_fluxes[flux_label][0]))
                    plt.plot(wavelength, flux_value*1E+3, 'o', color=color)

# Plot empirical fluxes for all subtypes
plot_empirical_fluxes(subtype_fluxes_and_colors_LDwarf , fluxes_LDwarf_empirical)
plot_empirical_fluxes(subtype_fluxes_and_colors_TDwarf, fluxes_TDwarf_empirical)
plot_empirical_fluxes(subtype_fluxes_and_colors_MDwarf, fluxes_MDwarf_empirical)

# Add labels and legend
plt.xlabel("Central Wavelength (Å)")
plt.ylabel("Flux (mJy)")
# plt.xlim(0, 50000)
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()
