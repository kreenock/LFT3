import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

df = pd.read_csv('sublunar_points.txt')
longs = df['Elong']
lats = df['Elat']
jds = df['JD']

interp = interp1d(jds, longs)

plt.plot(longs, lats, '.')
plt.show()