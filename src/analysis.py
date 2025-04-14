import numpy as np
from photutils.detection import DAOStarFinder

def calculate_flux(image_data, center, radius):
    """
    Функция, которая считает поток в заданной круговой апертуре в заданной точке.
    """
    y, x = np.ogrid[-float(center[1]):image_data.shape[0] - float(center[1]),
           -float(center[0]):image_data.shape[1] - float(center[0])]
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


def plot_graph_var(dataframe, filter1):
    """
    Функция, которая строит графики (время, звездная величина)
    """
    for index, row in dataframe.iterrows():
        num = row['Number']  # номер файла
        anchor1 = ast.literal_eval(row['Опорная1'])
        anchor2 = ast.literal_eval(row['Опорная2'])
        var1 = ast.literal_eval(row['Переменная1'])
        var2 = ast.literal_eval(row['Переменная2'])
        coords = [anchor1, anchor2, var1, var2]
        filename = '/home/vladislav/Education/3 курс/astro_prac/photometry_of_variable_stars/data/' + str(filter1) + '/' + str(filter1) + str(num) + '.fit'

        hdul = fits.open(filename)
        shot_data = hdul[0].data
        date_obs = hdul[0].header.get('DATE-OBS')
        obs_time = Time(date_obs, format='isot')
        our_times.append(obs_time.jd)
        hdul.close()
        flux = magn_finder(coords, shot_data, 7.80, 9.01, 5)
        our_magnit.append(flux)

    flux1, flux2 = zip(*our_magnit)

    magn1 = pd.Series(flux1)
    magn1.name = 'Flux1'
    magn2 = pd.Series(flux2)
    magn2.name = 'Flux2'
    dataframe = dataframe.join(magn1)
    dataframe = dataframe.join(magn2)
    dataframe = dataframe.sort_values(by='jd')

    max_flux_indices = dataframe['Flux2'].nlargest(2).index
    min_flux_indices = dataframe['Flux2'].nsmallest(3).index
    filtered_df = dataframe.drop(max_flux_indices)
    filtered_df = dataframe.drop(min_flux_indices)
    filtered_df = filtered_df[~((filtered_df['jd'] > 2.4535e6) & (filtered_df['jd'] < 2.4542e6))][['jd', 'Flux2']]
    jd = filtered_df['jd'].values
    flux2 = filtered_df['Flux2'].values

    spline2 = UnivariateSpline(jd, flux2, s=4)
    jd_smooth = np.linspace(jd.min(), jd.max(), 300)
    flux2_smooth = spline2(jd_smooth)

    plt.figure(figsize=(10, 6))
    plt.plot(jd, flux2, label='Flux 2 (not Original)', marker='o', linestyle='none')
    # plt.plot(df_h['jd'], df_h['Flux2'], label='Flux 2 (Original)', marker='o', linestyle='None')
    plt.plot(jd_smooth, flux2_smooth, label='Flux 2 (Smoothed)', linestyle='-')
    plt.xlabel('Time (Julian Date)')
    plt.ylabel('Flux')
    plt.title('Flux vs Time')
    plt.legend()
    plt.grid(True)
    plt.show()
