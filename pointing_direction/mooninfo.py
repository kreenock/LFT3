import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
from astroquery.jplhorizons import Horizons

df = pd.read_csv('sub2.csv')
longs = df['Elong']
lats = df['Elat']
jds = df['JD']

mobj = Horizons(id=301, location=500, epochs={'start': '2028-06-01T00:00:00', 'stop': '2028-06-30T23:59:59', 'step': '1h'})
meph = mobj.ephemerides()
obs_lon = meph['PDObsLon']
for i in range(0, len(obs_lon)):
    if obs_lon[i] > 180:
        obs_lon[i] -= 360
x_interp = interp1d(jds, longs)
y_interp = interp1d(jds, lats)

plt.plot(longs, lats, '.')
plt.plot(x_interp(np.arange(2461923.5, 2461953.0, 0.166)), y_interp(np.arange(2461923.5, 2461953.0, 0.166)), '.', color='r')
#plt.plot(obs_lon, meph['PDObsLat'], '.', color='r')
plt.xlabel('Selenographic longitude')
plt.ylabel('Selenographic latitude')
plt.title('Sublunar point June 2028; blue=grid method, red=interpolated')
plt.show()

# obj = Horizons(id=301, location=500, epochs={'start': '2028-06-01T0:00', 'stop': '2028-06-01T1:01', 'step':'30m'})
# eph = obj.ephemerides()
# print(eph['datetime_jd'])
