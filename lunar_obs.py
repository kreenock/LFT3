from astroquery.jplhorizons import Horizons
from pygdsm import GlobalSkyModel
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import healpy as hp
from geometry import angles

# see Danny Price https://github.com/telegraphic/pygdsm/tree/master
NSIDE = 512


class Galaxy:
    def __init__(self, freqs, fwhm, idisplay):
        self.gsm = GlobalSkyModel(include_cmb=True)  # This uses GSM08 (?)
        self.freqs = freqs
        self.fwhm = fwhm
        self.idisplay = idisplay

    def gen_map_cube(self):
        map_raw = self.gsm.generate(self.freqs)
        if self.fwhm is None:
            self.map_cube = map_raw
        else:
            self.map_cube = []
            for i in range(len(self.freqs)):
                self.map_cube.append(hp.sphtfunc.smoothing(map_raw[i], self.fwhm[i]))
        self.map_cube = np.array(self.map_cube)

    def gen_pointings(self):
        self.coord = Pointing(182.13737, -23.78930, '2028-01-01T00:00:00', '2028-12-31T23:59:59', '1h')
        self.locations = hp.pixelfunc.ang2pix(NSIDE, theta=self.coord.view.galactic.l.to_value(),
                                              phi=self.coord.view.galactic.b.to_value(), lonlat=True)

    def view(self, view_fullres=True, logged=True):
        plt.figure('Pointings')
        if view_fullres:
            self.gsm.view(self.idisplay, logged=True)
        else:
            if logged:
                hp.visufunc.mollview(np.log10(self.map_cube[self.idisplay]))
            else:
                hp.visufunc.mollview(self.map_cube[self.idisplay])


class Pointing:
    def __init__(self, telescope_long, telescope_lat, start_time, end_time, step):
        telescope = {'lon': telescope_long * u.deg, 'lat': telescope_lat * u.deg, 'elevation': 0 * u.km, 'body': 301}
        self.telescope = telescope
        self.start_time = start_time
        self.end_time = end_time
        self.step = step

        mobj = Horizons(id=301, location=500, epochs={'start': start_time, 'stop': end_time, 'step': step})
        meph = mobj.ephemerides()


        dec_difs = []
        self.get_difs(meph, dec_difs, telescope_long, telescope_lat)
        self.pointing_RAs = []
        self.pointing_DECs = []
        self.offset(self.pointing_RAs, self.pointing_DECs, dec_difs, meph['RA'], meph['DEC'])
        self.view = SkyCoord(ra=np.array(self.pointing_RAs)*u.degree, dec=np.array(self.pointing_DECs)*u.degree, frame='icrs')
        self.times = meph['datetime_jd']

        low_DECs, low_RAs = self.view_dif(self.pointing_RAs, self.pointing_DECs, -60)
        hi_DECs, hi_RAs = self.view_dif(self.pointing_RAs, self.pointing_DECs, 60)
        self.lowview = SkyCoord(ra=np.array(low_RAs) * u.degree, dec=np.array(low_DECs) * u.degree, frame='icrs')
        self.hiview = SkyCoord(ra=np.array(hi_RAs) * u.degree, dec=np.array(hi_DECs) * u.degree, frame='icrs')

    def get_difs(self, ephs, dec_list, long, lat):
        o_long, o_lat = self.get_origin()
        for i in range(0, len(ephs['DEC'])):
            dec_dif = angles(lat, long, o_lat[i], o_long[i])
            dec_list.append(dec_dif)

    def get_origin(self):
        mobj = Horizons(id=301, location=500,
                        epochs={'start': self.start_time, 'stop': self.end_time, 'step': self.step})
        meph = mobj.ephemerides()
        obs_lon = meph['PDObsLon']
        for i in range(0, len(obs_lon)):
            if obs_lon[i] > 180:
                obs_lon[i] -= 360
        return obs_lon, meph['PDObsLat']

    def offset(self, r_list, d_list, d_difs, ras, decs):
        for i in range(0, len(decs)):
            new_DEC = decs[i] + d_difs[i]
            new_RA = ras[i] #* np.cos(np.deg2rad(d_difs[i]))
            if new_DEC < -90:
                # dave, chatgpt, and i all had different opinions on this. chatgpt's method made the best looking plot
                # so this is chatgpt's solution for new_DEC < -90
                new_DEC = 90.0 + (new_DEC + 90)
                new_RA = ras[i] + 180
                if new_RA >= 360:
                    new_RA -= 360
            elif new_DEC > 90:
                new_DEC = 90 - (new_DEC - 90)
                new_RA = ras[i] + 180
                if new_RA >= 360:
                    new_RA -= 360
            d_list.append(new_DEC)
            r_list.append(new_RA)

    def view_dif(self, ras, decs, dif):
        new_DECs = np.array([])
        new_RAs = np.array([])
        for i in range(0, len(decs)):
            new_DEC = decs[i] + dif
            new_DECs = np.append(new_DECs, new_DEC)
            new_RAs = np.append(new_RAs, ras[i])

        low_flip = np.where(new_DECs < -90)
        # This is what I should do:
        # new_DECs[low_flip] = 90.0 - (new_DECs[low_flip] + 90)
        # new_RAs[low_flip] + 180

        # However, this flipped section ends up being completely inside the higher unviewed section, and plotting an
        # unviewed area within a bigger unviewed area is redundant. So I'm just going to remove it entirely.
        new_DECs = np.array([i for j, i in enumerate(new_DECs) if j not in low_flip[0]])
        new_RAs = np.array([i for j, i in enumerate(new_RAs) if j not in low_flip[0]])

        hi_flip = np.where(new_DECs > 90)
        new_DECs[hi_flip] = 90 - (new_DECs[hi_flip] - 90)
        new_RAs[hi_flip] += 180
        return new_DECs, new_RAs

