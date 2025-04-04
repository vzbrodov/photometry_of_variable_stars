from math import log10

from astropy.constants.codata2010 import alpha
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
from matplotlib.mlab import magnitude_spectrum
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std

def calculate_flux(image_data, center, radius):
    """
    Р¤СѓРЅРєС†РёСЏ СЃС‡РёС‚Р°РµС‚ Р·РІРµР·РґРЅСѓСЋ РІРµР»РёС‡РёРЅСѓ. РџСЂРёРЅРёРјР°РµС‚СЃСЏ РёР·РѕР±СЂР°Р¶РµРЅРёРµ, С†РµРЅС‚СЂ РіР°Р»Р°РєС‚РёРєРё, СЂР°РґРёСѓСЃ,
    РЅСѓР»СЊ-РїСѓРЅРєС‚ С€РєР°Р»С‹, РІСЂРµРјСЏ СЌРєСЃРїРѕР·РёС†РёРё, Р·РµРЅРёС‚РЅРѕРµ СЂР°СЃСЃС‚РѕСЏРЅРёСЏ, РєРѕСЌС„С„РёС†РёРµРЅС‚ РїРѕРіР»РѕС‰РµРЅРёСЏ РІ РґР°РЅРЅРѕР№ С†РІРµС‚РѕРІРѕР№ РїРѕР»РѕСЃРµ.
    Р’РѕР·РІСЂР°С‰Р°РµС‚ РѕРґРЅРѕ Р·РЅР°С‡РµРЅРёРµ -- Р·РІРµР·РґРЅСѓСЋ РІРµР»РёС‡РёРЅСѓ
    """

    y, x = np.ogrid[-center[1]:image_data.shape[0] - center[1], -center[0]:image_data.shape[1] - center[0]]

    mask = x ** 2 + y ** 2 <= radius ** 2

    flux = np.sum(image_data[mask])

    return flux

data_path_k = "C:/Users/Lenovo/Downloads/lab2/k/"
image_k = [
    data_path_k + "k60.fit"
]

hdul = fits.open(image_k[0])
data = hdul[0].data  # РџСЂРµРґРїРѕР»Р°РіР°РµРј, С‡С‚Рѕ РґР°РЅРЅС‹Рµ РІ РїРµСЂРІРѕРј HDU
hdul.close()


bkg_sigma = np.median(data) #mad_std(data)
print(bkg_sigma)
daofind = DAOStarFinder(fwhm = 5.0, threshold = 2000)
sources = daofind(data)


if sources:
    xy_coords = np.column_stack((sources['xcentroid'], sources['ycentroid']))

print(xy_coords)  # РњР°СЃСЃРёРІ shape (N, 2)


def magn_finder(coords,data,anchor1_mag,anchor2_mag,radius):
    #coords=(anch1, anch2, star1, star2)
    flux=[]

    anch1_flux = calculate_flux(data, coords[0], radius)
    anch2_flux = calculate_flux(data, coords[1], radius)
    star1_flux = calculate_flux(data, coords[2], radius)
    star2_flux = calculate_flux(data, coords[3], radius)

    flux.append( (anchor1_mag -2.5 * np.log10( star1_flux / anch1_flux ) + anchor2_mag -2.5 * np.log10( star1_flux / anch2_flux )) / 2 )
    flux.append( (anchor1_mag - 2.5 * np.log10(star2_flux / anch1_flux ) + anchor2_mag - 2.5 * np.log10(star2_flux / anch2_flux )) / 2 )

    return flux


print(magn_finder([[74.6674061 , 84.65724811],[176.18033724,  92.50422439],[176.18033724,  92.50422439],[75.50812106, 55.79575845]],data,7.42,6.93,9))