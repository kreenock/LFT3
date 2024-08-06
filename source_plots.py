from astroquery.jplhorizons import Horizons
import matplotlib.pyplot as plt
import numpy as np
from pointing_direction.get_pointing import Pointing
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd


class sources:
    def __init__(self, start_time, end_time, step):
        beam = Pointing(182.13737, -23.78930, start_time, end_time, step)
        sun = Horizons(id='sun', location=500,
                       epochs={'start': start_time, 'stop': end_time, 'step': step})
        jup = Horizons(id=599, location=500,
                       epochs={'start': start_time, 'stop': end_time, 'step': step})
        aRA = 219.90206
        aDEC = -60.83399
        stars = pd.read_fwf("closest_stars.txt", dtype=None)

        bRAs = beam.RAs
        bDECs = beam.DECs
        btimes = beam.times

        # for the other beams
        # for i in range(-5, 6, 1):
        #     new_lat = -23.78930 + (i * 10)
        #     if new_lat < -90:
        #         new_lat = -90 + ((-1 * new_lat) - 90)
        #     t = Pointing(182.13737, new_lat, '2028-01-01T00:00:00', '2028-12-31T23:59:59', '1h')
        #     bRAs.extend(t.RAs)
        #     bDECs.extend(t.DECs)
        #     btimes.extend(t.times)

        sRAs = sun.ephemerides()['RA']
        sDECs = sun.ephemerides()['DEC']
        stimes = sun.ephemerides()['datetime_jd']

        jRAs = jup.ephemerides()['RA']
        jDECs = jup.ephemerides()['DEC']
        jtimes = jup.ephemerides()['datetime_jd']

        spRAs = []
        spDECs = []
        sptimes = []

        for i in range(0, len(sRAs)):
            vis = False
            for j in np.where(btimes == stimes[i])[0]:
                if (bRAs[j] - 5) <= sRAs[i] <= (bRAs[j] + 5):
                    if (bDECs[j] - 60) <= sDECs[i] <= (bDECs[j] + 60):
                        vis = True
                        break
                if vis:
                    break
            if vis:
                spRAs.append(sRAs[i])
                spDECs.append(sDECs[i])
                sptimes.append(stimes[i])

        jpRAs = []
        jpDECs = []
        jptimes = []

        for i in range(0, len(jRAs)):
            vis = False
            for j in np.where(btimes == jtimes[i])[0]:
                if (bRAs[j] - 5) <= jRAs[i] <= (bRAs[j] + 5):
                    if (bDECs[j] - 60) <= sDECs[i] <= (bDECs[j] + 60):
                        vis = True
                        break
                if vis:
                    break
            if vis:
                jpRAs.append(jRAs[i])
                jpDECs.append(jDECs[i])
                jptimes.append(jtimes[i])

        apRAs = []
        apDECs = []
        aptimes = []

        for i in range(0, len(jtimes)):
            vis = False
            for j in np.where(btimes == jtimes[i])[0]:
                if (bRAs[j] - 5) <= aRA <= (bRAs[j] + 5):
                    if (bDECs[j] - 60) <= aDEC <= (bDECs[j] + 60):
                        vis = True
                        break
                if vis:
                    break
            if vis:
                apRAs.append(aRA)
                apDECs.append(aDEC)
                aptimes.append(jtimes[i])

        stpRAs = []
        stpDECs = []

        for i in range(0, len(stars['RA'])):
            vis = False
            for j in range(0, len(btimes)):
                if (bRAs[j] - 5) < stars['RA'][i] < (bRAs[j] + 5):
                    if (bDECs[j] - 60) < stars['DEC'][i] < (bDECs[j] + 60):
                        vis = True
                        break
                if vis:
                    break
            if vis:
                stpRAs.append(stars['RA'][i])
                stpDECs.append(stars['DEC'][i])

        self.sunview = SkyCoord(ra=np.array(spRAs)*u.degree, dec=np.array(spDECs)*u.degree, frame='icrs')
        self.jupview = SkyCoord(ra=np.array(jpRAs) * u.degree, dec=np.array(jpDECs) * u.degree, frame='icrs')
        self.aceview = SkyCoord(ra=np.array(apRAs) * u.degree, dec=np.array(apDECs) * u.degree, frame='icrs')
        self.staview = SkyCoord(ra=np.array(stpRAs) * u.degree, dec=np.array(stpDECs) * u.degree, frame='icrs')

        self.sun = [spRAs, spDECs, sptimes]
        self.jup = [jpRAs, jpDECs, jptimes]
        self.ace = [apRAs, apDECs, aptimes]
        self.sta = [spRAs, spDECs]

# print('Sun:', len(spRAs), spRAs)
# fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3)
#
# ax0.set_title('RA vs time')
# ax0.plot(btimes, bRAs, ',', color='b')
# ax0.plot(sptimes, spRAs, '.', color='y')
# ax0.plot(jptimes, jpRAs, '.', color='r')
# ax0.plot(aptimes, apRAs, '.', color='orange')
#
# ax1.set_title('Dec vs time')
# ax1.plot(btimes, bDECs, ',', color='b')
# ax1.plot(sptimes, spDECs, '.', color='y')
# ax1.plot(jptimes, jpDECs, '.', color='r')
# ax1.plot(aptimes, apDECs, '.', color='orange')
#
#
# ax2.set_title('RA vs Dec')
# ax2.plot(bRAs, bDECs, ',', color='b')
# ax2.plot(spRAs, spDECs, '.', color='y')
# ax2.plot(jpRAs, jpDECs, '.', color='r')
# ax2.plot(apRAs, apDECs, '.', color='orange')
#
# fig.suptitle('Visibility 2028')
# plt.show()
