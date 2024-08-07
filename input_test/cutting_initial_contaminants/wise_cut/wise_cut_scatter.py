import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Read data
data_qso_araa = ascii.read('QSO_wise_spec.dat')
data_qso_milliquas = ascii.read('QSO_milliquas_wise.dat')
data_qso_yang = ascii.read('QSO_yang_wise.dat')
data_bd_wise = ascii.read('BD_wise_SpT.dat', encoding='latin-1')
data_bd_temp_L = ascii.read("L_type_wise.txt")
data_bd_temp_M = ascii.read("M_type_wise.txt")
data_bd_temp_T = ascii.read("T_type_wise.txt")

# Function to convert column to float, replacing invalid entries with NaN
def convert_to_numeric(column):
    cleaned_column = [str(x).replace(',', '') for x in column]  # Remove commas
    return np.array(cleaned_column, dtype=float)

# Quasars ARAA
W1_qso_araa = data_qso_araa['w1mpro'] + 2.699
W2_qso_araa = data_qso_araa['w2mpro'] + 3.399

# Quasars Million Quasars Catalog v8
W1_qso_milliquas = data_qso_milliquas['w1mpro'] + 2.699
W2_qso_milliquas = data_qso_milliquas['w2mpro'] + 3.399

#Quasars Yang 2023
W1_qso_yang = data_qso_yang["w1mpro"] + 2.699
W2_qso_yang = data_qso_yang["w2mpro"] + 3.399

# Brown Dwarfs list (Bañados+2016)
W1_bd = data_bd_wise['w1mpro'] + 2.699
W2_bd = data_bd_wise['w2mpro'] + 3.399
Sp_T_bd = data_bd_wise['optical_type']

# Brown Dwarfs Empirical templates
Sp_T_bd_temp_L = data_bd_temp_L.columns[0]
Sp_T_bd_temp_M = data_bd_temp_M.columns[0]
Sp_T_bd_temp_T = data_bd_temp_T.columns[0]
W1_bd_temp_L = convert_to_numeric(data_bd_temp_L.columns[1]) + 2.699
W2_bd_temp_L = convert_to_numeric(data_bd_temp_L.columns[2]) + 3.399
W1_bd_temp_M = convert_to_numeric(data_bd_temp_M.columns[1]) + 2.699
W2_bd_temp_M = convert_to_numeric(data_bd_temp_M.columns[2]) + 3.399
W1_bd_temp_T = convert_to_numeric(data_bd_temp_T.columns[1]) + 2.699
W2_bd_temp_T = convert_to_numeric(data_bd_temp_T.columns[2]) + 3.399

# Filter the brown dwarfs by spectral type
M_type_mask = np.char.startswith(np.array(Sp_T_bd, dtype=str), 'M')
L_type_mask = np.char.startswith(np.array(Sp_T_bd, dtype=str), 'L')
T_type_mask = np.char.startswith(np.array(Sp_T_bd, dtype=str), 'T')
other_type_mask = ~M_type_mask & ~L_type_mask & ~T_type_mask

# Plotting
# Define the colors for each type
color_M = "purple"
color_L = "orange"
color_T = "black"

plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(W1_qso_araa, W1_qso_araa - W2_qso_araa, label="Quasars Fan+2022", s=10, alpha=0.4, marker="x", color="green")
plt.scatter(W1_qso_milliquas, W1_qso_milliquas - W2_qso_milliquas, label="Quasars Flesch+2023", s=10, alpha=0.2, marker="x", color="blue")
plt.scatter(W1_qso_yang, W1_qso_yang - W2_qso_yang, label="Quasars Yang+2023", s=10, alpha=0.4, marker="x", color="red")
plt.scatter(W1_bd[other_type_mask], W1_bd[other_type_mask] - W2_bd[other_type_mask], label="Brown Dwarfs without SpT", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(W1_bd[M_type_mask], W1_bd[M_type_mask] - W2_bd[M_type_mask], label="M-type Brown Dwarf", s=8, alpha=0.2, marker="o", color=color_M)
plt.scatter(W1_bd[L_type_mask], W1_bd[L_type_mask] - W2_bd[L_type_mask], label="L-type Brown Dwarf", s=8, alpha=0.2, marker="o", color=color_L)
plt.scatter(W1_bd[T_type_mask], W1_bd[T_type_mask] - W2_bd[T_type_mask], label="T-type Brown Dwarf", s=8, alpha=0.8, zorder=1, marker="o", color=color_T)
plt.scatter(W1_bd_temp_T, W1_bd_temp_T - W2_bd_temp_T, label="T-type Brown Dwarf Template", s=100, marker="*", color=color_T)
plt.scatter(W1_bd_temp_M, W1_bd_temp_M - W2_bd_temp_M, label="M-type Brown Dwarf Template", s=100, marker="*", color=color_M)
plt.scatter(W1_bd_temp_L, W1_bd_temp_L - W2_bd_temp_L, label="L-type Brown Dwarf Template", s=100, marker="*", color=color_L)
plt.xlabel("W1")
plt.ylabel("W1-W2")
plt.axhline(y=-0.6, color='gray', linestyle='--', label='WISE color cut')
plt.axhline(y=0.6, color='gray', linestyle='--')
plt.legend(loc="upper left")
plt.ylim(-1, 2.5)
# plt.show()
plt.savefig("QSO_BD_wise_colors.png")