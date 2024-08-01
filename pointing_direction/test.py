from astroquery.jplhorizons import Horizons
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from get_pointing import Pointing
import pandas as pd

stars = pd.read_fwf("closest_stars.txt", dtype=None)

