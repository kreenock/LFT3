from astroquery.jplhorizons import Horizons
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np


def get_vectors(target, observer, t):
    obj = Horizons(id=target, location=observer, epochs=t)
    return obj.vectors()


def get_origin(time, start_lat, stop_lat, start_long, stop_long, step_lat, step_long):
    o_list = []
    lat_list = []
    long_list = []
    for i in np.arange(start_long, (stop_long + 1), step_long):
        for j in np.arange(start_lat, (stop_lat + 1), step_lat):
            origin_test = {'lon': i * u.deg, 'lat': j * u.deg, 'elevation': 0 * u.km, 'body': 301}
            vecs = get_vectors(origin_test, 500, time)
            distance_km = (vecs['x'] ** 2 + vecs['y'] ** 2 + vecs['z'] ** 2)
            o_list.append(distance_km)
            lat_list.append(j)
            long_list.append(i)
    return lat_list[o_list.index(min(o_list))], long_list[o_list.index(min(o_list))]


o_lats = []
o_longs = []

timer = 2460310.5
for i in range(0, 744):
    if i == 0:
        o_lat, o_long = get_origin(str(timer), -10, 10, -10, 10, 1, 1)
    else:
        o_lat, o_long = get_origin(str(timer), (nola - 2), (nola + 2), (nolo - 2), (nolo + 2), 0.5, 0.5)
    o_lats.append(o_lat)
    o_longs.append(o_long)
    nola = o_lat
    nolo = o_long
    timer += 0.04166666667

plt.figure()
plt.scatter(o_lats, o_longs)
plt.title("Sublunar point every hour June 2028")
plt.xlabel("Longitude (degrees)")
plt.ylabel("Latitude (degrees")
print(o_lats, o_longs)
plt.show()
