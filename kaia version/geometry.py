import numpy as np


def angles(tele_lat, tele_long, origin_lat, origin_long):
    ctele_lat = 90 - tele_lat
    tele_x = np.sin(np.deg2rad(ctele_lat)) * np.cos(np.deg2rad(tele_long))
    tele_y = np.sin(np.deg2rad(ctele_lat)) * np.sin(np.deg2rad(tele_long))
    tele_z = np.cos(np.deg2rad(ctele_lat))
    tele_mag = (tele_x ** 2 + tele_y ** 2 + tele_z ** 2) ** 0.5

    anti_lat = -1 * origin_lat

    anti_long = origin_long + 180
    if anti_long >= 360:
        anti_long -= 360

    canti_lat = 90 - anti_lat

    anti_x = np.sin(np.deg2rad(canti_lat)) * np.cos(np.deg2rad(anti_long))
    anti_y = np.sin(np.deg2rad(canti_lat)) * np.sin(np.deg2rad(anti_long))
    anti_z = np.cos(np.deg2rad(canti_lat))
    anti_mag = (anti_x ** 2 + anti_y ** 2 + anti_z ** 2) ** 0.5

    dot = (tele_x * anti_x) + (tele_y * anti_y) + (tele_z * anti_z)
    theta = (np.arccos(dot / (tele_mag * anti_mag)) * 180 / np.pi)

    if tele_lat < anti_lat:
        theta = theta * -1

    return theta
