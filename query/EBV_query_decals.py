import numpy as np
import dustmaps.sfd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from dustmaps.sfd import SFDQuery
from dustmaps.config import config



# Configure dustmaps data directory
config['data_dir'] = '/Users/dathev/anaconda3/lib/python3.10/site-packages/dustmaps'
dustmaps.sfd.fetch()

# Load the SFD dust map
sfd = SFDQuery()

# Read the .fits table
input_filename = 'output_catalogs/results_catalogs/final_TM.fits'  # Replace with your input filename
output_filename = 'final_TM_ebv'  # Replace with your output filename
with fits.open(input_filename) as hdul:
    data = hdul[1].data

# Prepare the new column for E(B-V)
ebv_values = np.zeros(len(data))

# Query E(B-V) for each source
for i, row in enumerate(data):
    coords = SkyCoord(ra=row['ra'], dec=row['dec'], unit='deg', frame='icrs')
    ebv_values[i] = sfd(coords)

# Add the new column to the table
table = Table(data)
table['E(B-V)'] = ebv_values

# Write the new table to a .fits file
table.write(output_filename, format='fits', overwrite=True)

print(f"Output saved to {output_filename}")