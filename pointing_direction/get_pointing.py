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
        print(meph['datetime_jd'][0])

        dec_difs = []
        self.get_difs(meph, dec_difs, telescope_long, telescope_lat)
        pointing_RAs = []
        pointing_DECs = []
        self.offset(pointing_RAs, pointing_DECs, dec_difs, meph['RA'], meph['DEC'])
        self.RAs = pointing_RAs
        self.DECs = pointing_DECs

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
            new_RA = ras[i] * np.cos(np.deg2rad(d_difs[i]))
            if new_DEC < -90:
                new_DEC = -90 + ((-1 * new_DEC) - 90)
            elif new_DEC > 90:
                new_DEC = 90 - (new_DEC - 90)
            d_list.append(new_DEC)
            r_list.append(new_RA)


if __name__ == "__main__":
    tele = Pointing(182.13737, -23.78930, '2028-06-01T0:00:00', '2028-06-30T23:59:59', '1h')
    RAs = tele.RAs
    DECs = tele.DECs
    print(DECs, RAs)
    plt.plot(RAs, DECs, '.')
    plt.title('Pointing RA vs Dec June 2028')
    plt.xlabel('Dec')
    plt.ylabel('RA')
#    plt.savefig("pointing_plot.jpg")
    plt.show()