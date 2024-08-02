from astroquery.jplhorizons import Horizons
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import numpy as np
import matplotlib.gridspec as gridspec

sun = Horizons(id='sun', location='500', epochs={'start': '2027-01-01T00:00', 'stop': '2027-12-31T23:59', 'step': '24h'})
seph = sun.ephemerides()
jupiter = Horizons(id='599', location='500', epochs={'start': '2027-01-01T00:00', 'stop': '2027-12-31T23:59', 'step': '24h'})
jeph = jupiter.ephemerides()

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
ax.scatter(seph['RA']/30, seph['DEC'], c='y')
ax.scatter(jeph['RA']/30, jeph['DEC'], c='r')
plt.title("theta=RA, r=Dec")
plt.show()