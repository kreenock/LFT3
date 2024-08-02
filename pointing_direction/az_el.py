from astroquery.jplhorizons import Horizons
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

telescope = {'lon': 182.13737, 'lat': -23.78930, 'elevation': 0, 'body': 301}
start_time = '2028-06-01T00:00'
end_time = '2028-06-31T23:59'

sobj = Horizons(id='sun', location=301, epochs={'start': start_time, 'stop': end_time, 'step': '1h'})
seph = sobj.ephemerides()

jobj = Horizons(id='599', location=301, epochs={'start': start_time, 'stop': end_time, 'step': '1h'})
jeph = jobj.ephemerides()

ss = np.where(np.array(seph['EL']) >= 30)[0]
sazs = []
sels = []
sts = []
for i in ss:
    sazs.append(seph['RA'][i])
    sels.append(seph['DEC'][i])
    sts.append(seph['datetime_jd'][i])

jj = np.where(np.array(jeph['EL']) >= 30)[0]
jazs = []
jels = []
jts = []
for j in jj:
    jazs.append(jeph['RA'][j])
    jels.append(jeph['DEC'][j])
    jts.append(jeph['datetime_jd'][j])


fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.plot(sts, sazs, '.', color='y')
ax2.plot(sts, sels, '.', color='y')
ax3.plot(sazs, sels, '.', color='y')
ax1.plot(jts, jazs, '.', color='b')
ax2.plot(jts, jels, '.', color='b')
ax3.plot(jazs, jels, '.', color='b')
ax1.set_title('Sun/Jupiter RA vs time')
ax2.set_title('Sun/Jupiter Dec vs time')
ax3.set_title('Sun/Jupiter RA vs Dec')
ax1.set_xlabel('Time (jd)')
ax1.set_ylabel('Azimuth (degrees)')
ax2.set_xlabel('Time (jd)')
ax2.set_ylabel('Elevation (degrees)')
ax3.set_xlabel('Azimuth (degrees)')
ax3.set_ylabel('Elevation (degrees)')
fig.suptitle('June 2028')


print(max(seph['datetime_jd']), max(jeph['datetime_jd']))
plt.show()
