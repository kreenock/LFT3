from astroquery.jplhorizons import Horizons
from astropy import units as u
import matplotlib.pyplot as plt


def pos_plot(target, observer, epoch, color):
    eph = Horizons(id=target, location=observer, epochs=epoch).ephemerides()
    plt.plot(eph['RA'], eph['DEC'], c=color)
#    plt.plot(eph['RA'][0], eph['DEC'][0], 'go')
#    plt.plot(eph['RA'][len(eph['RA']) - 1], eph['DEC'][len(eph['DEC']) - 1], 'ro')

    for i in range(0, len(eph)):
        if eph['datetime_str'][i][9: 11] == '01':
            plt.plot(eph['RA'][i], eph['DEC'][i], 'go')

    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")


start_time = '2028-01-01T00:00'
end_time = '2028-01-31T23:59'
times = {'start': start_time, 'stop': end_time, 'step': '12h'}
telescope = {'lon': 182.13737 * u.deg, 'lat': -23.78930 * u.deg, 'elevation': 0 * u.km, 'body': 301}

pos_plot(301, 500, times, 'y')
pos_plot(telescope, 500, times, 'b')

plt.show()
