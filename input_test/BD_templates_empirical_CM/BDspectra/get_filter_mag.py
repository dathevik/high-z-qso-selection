import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from synphot import SourceSpectrum, SpectralElement, Observation
from synphot.models import Empirical1D
from astropy.io import ascii
from astropy import constants as const
from astropy import units as u

# ===========================================================================

EXAMPLES = """
get_filter.py -i input_text_file -f folder name
"""


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='''
    Get a spectrum and convolve to filter to get magnitude. Use the Synphot package
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-i', '--input', required=True, type=str,
                        help='Path to the folder containing the spectrum files. The files should be ascii files with columns: wave (in micron) flux (in erg/cm2/s/micron).')

    parser.add_argument('-f', '--filter', required=True, type=str,
                        help='Filter name. Needs to have response in subdirectory filter/filter_XX.dat')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    # ===========================================================================
    # --- path where filter curves are.
    filter_name = args.filter + ".dat"
    bp_filt = SpectralElement.from_file(filter_name)

    # Filters in Angstrom
    bp_filt_wavelengths_ang= bp_filt.waveset.to(u.angstrom)

    # Create a new SpectralElement with wavelengths in microns
    bp_filt_ang = SpectralElement(Empirical1D, points=bp_filt_wavelengths_ang,
                                  lookup_table=bp_filt(bp_filt.waveset))

    # ===========================================================================
    # --- Process each file in the input directory
    input_folder = args.input
    output_filename = 'output_magnitudes.txt'
    cent_wav = 7827.66 * u.angstrom
    with open(output_filename, 'w') as output_file:
        output_file.write('Type\tMagnitude (VEGA)\tFlux (mJy)\n')

        for filename in os.listdir(input_folder):
            if filename.endswith('.txt'):
                filepath = os.path.join(input_folder, filename)
                temp_data = ascii.read(filepath)
                temp_wave = (temp_data.columns[0] * 10000) * u.angstrom
                temp_flux = (temp_data.columns[1] / 10000) * u.erg / (u.cm ** 2 * u.s * u.angstrom)

                # Ensure the source spectrum wavelengths are in microns
                temp_wave_ang = temp_wave.to(u.angstrom)

                # Check if the wavelength ranges overlap
                if (temp_wave_ang.max() < bp_filt_ang.waveset.min()) or (temp_wave_ang.min() > bp_filt_ang.waveset.max()):
                    print(f"Warning: The wavelength range of the file {filename} does not overlap with the filter.")
                    output_file.write(f'{filename[:-4]}\t\t\n')
                    continue

                # ---- Create spectrum object
                sp_temp = SourceSpectrum(Empirical1D, points=temp_wave_ang, lookup_table=temp_flux, keep_neg=True)

                # ---- Obtain Fluxes/Magnitudes
                obs_temp_bp = Observation(sp_temp, bp_filt_ang, binset=bp_filt_ang.waveset, force='extrap')
                flux = obs_temp_bp.effstim(u.erg / (u.cm ** 2 * u.s * u.angstrom)).value

                # Get Vega magnitude using synphot's built-in Vega spectrum
                vega_spectrum = SourceSpectrum.from_vega()
                vega_obs = Observation(vega_spectrum, bp_filt_ang)
                zero_point_flux = vega_obs.effstim(u.erg / (u.cm ** 2 * u.s * u.angstrom)).value

                vega_mag = -2.5 * np.log10(flux / zero_point_flux)
                flux_cgs = flux * 0.013616591668795114 * 1E-8
                flux_mjy = (3.33564095E+04 * flux_cgs * cent_wav ** 2) * 1E-3
                # Write results to the output file
                output_file.write(f'{filename[:-4]}\t{vega_mag:.4f}\t{flux_mjy:.4e}\n')

    print("===============================================================")
    print(f"Output file created: {output_filename}")
