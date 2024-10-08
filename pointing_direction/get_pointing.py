import numpy as np
from astroquery.jplhorizons import Horizons
from astropy import units as u
from geometry import angles
import matplotlib.pyplot as plt


class Pointing:
    def __init__(self, telescope_long, telescope_lat, start_time, end_time, step):
        telescope = {'lon': telescope_long * u.deg, 'lat': telescope_lat * u.deg, 'elevation': 0 * u.km, 'body': 301}
        self.telescope = telescope
        self.start_time = start_time
        self.end_time = end_time

        mobj = Horizons(id=301, location=500, epochs={'start': start_time, 'stop': end_time, 'step': step})
        meph = mobj.ephemerides()


        dec_difs = []
        self.get_difs(meph, dec_difs, telescope_long, telescope_lat)
        pointing_RAs = []
        pointing_DECs = []
        self.offset(pointing_RAs, pointing_DECs, dec_difs, meph['RA'], meph['DEC'])
        self.RAs = pointing_RAs
        self.DECs = pointing_DECs
        self.times = meph['datetime_jd']

    def get_difs(self, ephs, dec_list, long, lat):
        for i in range(0, len(ephs['DEC'])):
            #           o_lat, o_long = self.get_origin(ephs['datetime_jd'][i])
            o_lat = 0
            o_long = 0
            dec_dif = angles(lat, long, o_lat, o_long)
            dec_list.append(dec_dif)

    def get_origin(self, time):
        o_list = []
        lat_list = []
        long_list = []
        for i in np.arange(-9, 9, 0.5):
            for j in np.arange(-9, 9, 0.5):
                origin_test = {'lon': i, 'lat': j, 'elevation': 0, 'body': 301}
                obj = Horizons(id=origin_test, location=399, epochs=time)
                vecs = obj.vectors()
                distance_km = (vecs['x'] ** 2 + vecs['y'] ** 2 + vecs['z'] ** 2) ** 0.5
                o_list.append(distance_km)
                lat_list.append(j)
                long_list.append(i)
        return lat_list[o_list.index(min(o_list))], long_list[o_list.index(min(o_list))]

    def offset(self, r_list, d_list, d_difs, ras, decs):
        for i in range(0, len(decs)):
            new_DEC = decs[i] + d_difs[i]
            new_RA = ras[i] #* np.cos(np.deg2rad(d_difs[i]))
            if new_DEC < -90:
                new_DEC = 90.0 + (new_DEC + 90)
            elif new_DEC > 90:
                new_DEC = 90 - (new_DEC - 90)
            d_list.append(new_DEC)
            r_list.append(new_RA)


if __name__ == "__main__":
    tele = Pointing(182.13737, -23.78930, '2028-01-01T00:00:00', '2028-12-31T23:59:59', '1h')
    ras = tele.RAs
    decs = tele.DECs
    times = tele.times
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    ax1.plot(times, ras, '.')
    ax2.plot(times, decs, '.')
    ax3.plot(decs, ras, '.')
    ax1.set_title('Pointing RA vs time')
    ax2.set_title('Pointing Dec vs time')
    ax3.set_title('Pointing RA vs Dec')
    ax1.set_xlabel('Time (jd)')
    ax1.set_ylabel('RA (degrees)')
    ax2.set_xlabel('Time (jd)')
    ax2.set_ylabel('Dec (degrees)')
    ax3.set_xlabel('Dec (degrees)')
    ax3.set_ylabel('RA (degrees)')
    fig.suptitle('2028')
    plt.show()