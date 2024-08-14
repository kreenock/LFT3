import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.gaia import Gaia
import time
import sys
import argparse

#  --- CLI Arguments ---
parser = argparse.ArgumentParser(description='Gaia DR3 Cone Search')
parser.add_argument('-b', '--beamsize', default=60, type=float,
                    help='Beam size to query. Default is 2.59/2 degrees. Which is the FWHM of a LOFAR international '
                         'station at 150 MHz.')
parser.add_argument('-i', '--input', type=str, required=True,
                    help='Input file containing the RA and DEC of the zenith pointings. In the form of RA, '
                         'DEC. Delimited by a comma.')
args = parser.parse_args()

log_file = open('all-gaia-output.log', 'w')


class Tee(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()  # Flush the buffer to ensure immediate writing

    def flush(self):
        for f in self.files:
            f.flush()


# Create an instance of the Tee class to redirect the output
tee = Tee(sys.stdout, log_file)
sys.stdout = tee

# --- Setup Parameters ---
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Select Data Release 3
Gaia.ROW_LIMIT = -1  # No limit on number of rows returned

# --- Loading LOFAR beam pointings ---
OBS_df = pd.read_csv(args.input, delimiter=',', skiprows=1)
print(OBS_df.head())
pointings_ra = OBS_df['ra']
pointings_dec = OBS_df['dec']
pointings_vec = np.vstack((pointings_ra, pointings_dec)).T
pointing_coords = SkyCoord(ra=pointings_vec[:, 0] * u.degree, dec=pointings_vec[:, 1] * u.degree)
print('Number of pointings: ', len(pointings_vec))
beam_radius = args.beamsize

data_frames = []
total_count = 0
actual_count = 0
count_2σ = 0
count_3σ = 0

# Loop through the pointing coordinates and retrieve data from Gaia
for i in (range(0, len(pointing_coords))):
    start = time.time()
    # --- Query Gaia ---
    j = Gaia.cone_search_async(pointing_coords[i], radius=u.Quantity(beam_radius, u.deg))
    r = j.get_results()
    print('\n---\nTotal number of targets found in beam %s: %s' % (i, len(r)))
    print('RA (deg): ' '{:.3f}'.format(pointing_coords[i].ra.degree),
          '\nDEC (deg): ' '{:.3f}'.format(pointing_coords[i].dec.degree))
    total_count += len(r)
    df = r.to_pandas()

    # --- Seperation between pointing and target ---
    sep = pointing_coords[i].separation(SkyCoord(ra=df['ra'] * u.degree, dec=df['dec'] * u.degree)).degree
    coord_error = np.sqrt(df['ra_error'] ** 2 + df['dec_error'] ** 2)
    extension = sep + coord_error
    extension_2σ = sep + 2 * coord_error
    extension_3σ = sep + 3 * coord_error
    filter_mask_1σ = extension < beam_radius
    filter_mask_2σ = extension_2σ < beam_radius
    filter_mask_3σ = extension_3σ < beam_radius

    # --- apply filter mask ---
    print('Number of targets in beam within 1σ: ', len(df[filter_mask_1σ]))
    print('Number of targets in beam within 2σ: ', len(df[filter_mask_2σ]))
    print('Number of targets in beam within 3σ: ', len(df[filter_mask_3σ]))
    count_2σ += len(df[filter_mask_2σ])
    count_3σ += len(df[filter_mask_3σ])

    df = df[filter_mask_1σ]
    actual_count += len(df)
    end = time.time()
    print('Time taken to query Gaia for this beam: %s secs.' % "{:.1f}".format(end - start))

# Concatenate the individual data frames into a single total data frame
total_dataframe = pd.concat(data_frames, ignore_index=True)
total_dataframe.to_csv('GDR3_beam_targets_r%s_%s.csv' % (beam_radius, time.strftime("%Y%m%d")))

print('\n---\nTotal number of targets found in beam: ', total_count)
print('Number of targets found in beam (1σ): ', len(total_dataframe))
print('Number of targets found in beam (2σ): ', count_2σ)
print('Number of targets found in beam (3σ): ', count_3σ)
