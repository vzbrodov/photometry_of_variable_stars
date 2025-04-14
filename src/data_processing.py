import os
from photutils.detection import DAOStarFinder
import re
from astropy.io import fits
import numpy as np

def get_coords(data_path, filter1, my_threshold, my_fwhm):
    data_path = os.path.join(data_path, filter1) + '/'
    iterator_images = sorted(
        [os.path.join(data_path, f) for f in os.listdir(data_path) if f.endswith('.fit')],
        key=lambda x: int(re.search(rf'{filter1}(\d+)\.fit', x).group(1))
    )
    approx_coords = []

    for shot in iterator_images:
        print(f"Processing {shot}")
        hdul = fits.open(shot)
        shot_data = hdul[0].data
        hdul.close()

        bkg_sigma = np.median(shot_data)

        daofind = DAOStarFinder(fwhm=my_fwhm, threshold= my_threshold )
        sources = daofind(shot_data - bkg_sigma)  # Вычитаем фон

        if sources:
            xy_coords = np.column_stack((sources['xcentroid'], sources['ycentroid']))
            approx_coords.append(xy_coords)
            approx_coords.append(['==========','=========='])

    filename = 'OUTPUT' + filter1 + '.txt'
    with open(filename, 'w') as file:
        for coords in approx_coords:
            for coord in coords:
                file.write(f"{coord[0]} {coord[1]}\n")


def get_times(data_path, filter1):
    """
    Функция, которая получает времена наблюдения со всех снимков.
    """
    data_path = os.path.join(data_path, filter1) + '/'
    iterator_images = sorted(
        [os.path.join(data_path, f) for f in os.listdir(data_path) if f.endswith('.fit')],
        key=lambda x: int(re.search(rf'{filter1}(\d+)\.fit', x).group(1))
    )

    for file in iterator_images:
        with fits.open(file) as hdul:
            header = hdul[0].header

            date_obs = header.get('DATE-OBS')  # Дата наблюдения

            obs_time = Time(date_obs, format='isot')

            data.append({
                'observation_time': obs_time.isot,
                'jd': obs_time.jd  # Юлианская дата для расчетов
            })
