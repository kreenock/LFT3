from astroquery.jplhorizons import Horizons
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from copy import copy
import lunar_sys
import lunar_obs
import healpy as hp
from astropy.coordinates import SkyCoord
from source_plots import sources

NSIDE = 512
plt.style.use('ggplot')

color_palette = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0),
    (1.0, 0.4980392156862745, 0.054901960784313725, 1.0),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313, 1.0),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0),
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0),
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354, 1.0)
]


class Observe:
    def __init__(self, N=50, deck_diameter=3.0, element_low=400.0, array_low=None, start=250.0, stop=750, step=10):
        self.LB = [1, 50]
        self.FM = [70, 120]
        self.MA = [start, stop]
        self.HB = [start * 2, start * 6]
        self.system = lunar_sys.System(N=N, deck_diameter=deck_diameter, element_low=element_low, array_low=array_low,
                                       start=start, stop=stop, step=step)
        self.system.gen_Trcvr()
        self.system.gen_FWHM()
        self.system.check()
        self.Tsys = None
        self.galaxy = None

    def get_galaxy(self, pointing='moon_ptg'):
        self.galaxy = lunar_obs.Galaxy(freqs=self.system.freqs, fwhm=self.system.fwhm, idisplay=self.system.idisplay)
        self.galaxy.pointing_type = pointing
        self.galaxy.gen_map_cube()
        if pointing == 'moon_ptg':
            self.galaxy.gen_pointings()
        else:
            self.galaxy.gen_locs(112000)  # randomish
        self.galaxy.view()

    def get_sky_Tsys(self, pointing='moon_ptg'):
        if self.galaxy is None:
            self.get_galaxy(pointing=pointing)
        self.Gal = []
        self.Tsys = []
        for i in range(len(self.galaxy.locations)):
            self.Gal.append(self.galaxy.map_cube[:, self.galaxy.locations[i]])
            self.Tsys.append(self.galaxy.map_cube[:, self.galaxy.locations[i]] + self.system.Trcvr)

    def run(self):
        self.get_sky_Tsys()
        x = sources('2028-01-01T00:00:00', '2028-12-31T23:59:59', '1h')
        hp.visufunc.projplot(self.galaxy.coord.view.galactic.l.to_value(), self.galaxy.coord.view.galactic.b.to_value(),
                             'k', linestyle=':', lonlat=True)
        hp.visufunc.projplot(self.galaxy.coord.lowview.galactic.l.to_value(),
                             self.galaxy.coord.lowview.galactic.b.to_value(), 'k', lonlat=True)
        hp.visufunc.projplot(self.galaxy.coord.hiview.galactic.l.to_value(),
                             self.galaxy.coord.hiview.galactic.b.to_value(), 'k', lonlat=True)
        hp.visufunc.projplot(x.staview.galactic.l.to_value(), x.staview.galactic.b.to_value(), 'wx',
                             lonlat=True)
        hp.visufunc.projplot(x.sunview.galactic.l.to_value(), x.sunview.galactic.b.to_value(), 'yo',
                             lonlat=True)
        hp.visufunc.projplot(x.jupview.galactic.l.to_value(), x.jupview.galactic.b.to_value(), 'ro',
                             lonlat=True)
        hp.visufunc.projplot(x.marview.galactic.l.to_value(), x.marview.galactic.b.to_value(), 'co',
                             lonlat=True)
        hp.visufunc.projplot(x.aceview.galactic.l.to_value(), x.aceview.galactic.b.to_value(), 'go',
                             lonlat=True)