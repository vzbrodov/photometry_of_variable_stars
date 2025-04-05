import numpy as np
from photutils.detection import DAOStarFinder

def calculate_flux(image_data, center, radius):
    """
    Функция, которая считает поток в заданной круговой апертуре в заданной точке.
    """
    y, x = np.ogrid[-center[1]:image_data.shape[0] - center[1], -center[0]:image_data.shape[1] - center[0]]
    mask = x ** 2 + y ** 2 <= radius ** 2
    flux = np.sum(image_data[mask])
    return flux


def magn_finder(coords,data,anchor1_mag,anchor2_mag,radius):
    """
    Функция, которая считает звездные величины двух переменных звезд относительно
    двух опорных звезд, звездные величины которых известны заранее.
    На вход подаются координаты 4х звезд: (anch1, anch2, star1, star2)
    """
    #coords=(anch1, anch2, star1, star2)
    our_magnitudes = []

    anch1_flux = calculate_flux(data, coords[0], radius)
    anch2_flux = calculate_flux(data, coords[1], radius)
    star1_flux = calculate_flux(data, coords[2], radius)
    star2_flux = calculate_flux(data, coords[3], radius)

    our_magnitudes.append((anchor1_mag -2.5 * np.log10( star1_flux / anch1_flux ) + anchor2_mag -2.5 * np.log10( star1_flux / anch2_flux )) / 2 )
    our_magnitudes.append((anchor1_mag - 2.5 * np.log10(star2_flux / anch1_flux ) + anchor2_mag - 2.5 * np.log10(star2_flux / anch2_flux )) / 2 )

    return our_magnitudes
