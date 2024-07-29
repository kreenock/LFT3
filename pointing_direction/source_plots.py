import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from get_pointing import Pointing

tele = Pointing(182.13737, -23.78930, '2028-01-01T00:00:00', '2028-12-31T23:59:59', '1h')
sun = Horizons(id='sun', location=500, epochs={'start': '2028-01-01T00:00:00', 'stop': '2028-12-31T23:59:59', 'step': '1h'})
jup = Horizons(id=599, location=500, epochs={'start': '2028-01-01T00:00:00', 'stop': '2028-12-31T23:59:59', 'step': '1h'})

tRAs = tele.RAs
tDECs = tele.DECs
times = tele.times

sRAs = sun.ephemerides()['RA']
sDECs = sun.ephemerides()['DEC']

jRAs = jup.ephemerides()['RA']
jDECs = jup.ephemerides()['DEC']

spRAs = []
spDECs = []
sptimes = []

for i in range(0, len(sRAs)):
    if (tRAs[i] - 10) <= sRAs[i] <= (tRAs[i] + 10):
        if (tDECs[i] - 60) <= sDECs[i] <= (tDECs[i] + 60):
            spRAs.append(sRAs[i])
            spDECs.append(sDECs[i])
            sptimes.append(times[i])

jpRAs = []
jpDECs = []
jptimes = []

for i in range(0, len(jRAs)):
    if (tRAs[i] - 10) <= jRAs[i] <= (tRAs[i] + 10):
        if (tDECs[i] - 60) <= jDECs[i] <= (tDECs[i] + 60):
            jpRAs.append(jRAs[i])
            jpDECs.append(jDECs[i])
            jptimes.append(times[i])

fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3)

ax0.set_title('RA vs time')
ax0.plot(times, tRAs, '.', color='b')
ax0.plot(sptimes, spRAs, '.', color='y')
ax0.plot(jptimes, jpRAs, '.', color='r')

ax1.set_title('Dec vs time')
ax1.plot(times, tDECs, '.', color='b')
ax1.plot(sptimes, spDECs, '.', color='y')
ax1.plot(jptimes, jpDECs, '.', color='r')


ax2.set_title('RA vs Dec')
ax2.plot(tRAs, tDECs, '.', color='b')
ax2.plot(spRAs, spDECs, '.', color='y')
ax2.plot(jpRAs, jpDECs, '.', color='r')

fig.suptitle('Visibility 2028')
plt.show()
