from astroquery.jplhorizons import Horizons
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import math

telescope = {'lon': 182.13737 * u.deg, 'lat': -23.78930 * u.deg, 'elevation': 0 * u.km, 'body': '301'}
times = {'start': '2028-01-07T00:00', 'stop': '2028-1-21T23:59', 'step': '24h'}


def make_plot(target, observer, epochs, color):
    obj = Horizons(id=target, location=observer, epochs=epochs)
    eph = obj.ephemerides()

    plt.plot(eph['RA'], eph['DEC'], c=color)
    plt.plot(eph['RA'][0], eph['DEC'][0], 'go')
    plt.plot(eph['RA'][len(eph['RA']) - 1], eph['DEC'][len(eph['DEC']) - 1], 'ro')


make_plot(telescope, 301, times, 'y')
make_plot(301, 399, times, 'b')
plt.show()