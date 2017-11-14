
from wm_nctools import save_multiplier
import numpy as np


from ProcessMultipliers import processMultipliers as pM

indices = {
    0: {'dir': 'n', 'min': 0., 'max': 22.5, 'fill': 0},
    1: {'dir': 'ne', 'min': 22.5, 'max': 67.5, 'fill': 1},
    2: {'dir': 'e', 'min': 67.5, 'max': 112.5, 'fill': 2},
    3: {'dir': 'se', 'min': 112.5, 'max': 157.5, 'fill': 3},
    4: {'dir': 's', 'min': 157.5, 'max': 202.5, 'fill': 4},
    5: {'dir': 'sw', 'min': 202.5, 'max': 247.5, 'fill': 5},
    6: {'dir': 'w', 'min': 247.5, 'max': 292.5, 'fill': 6},
    7: {'dir': 'nw', 'min': 292.5, 'max': 337.5, 'fill': 7},
    8: {'dir': 'n', 'min': 337.5, 'max': 360., 'fill': 0},
    9: {'dir': 'max'}
}

def example():
    multiplier_name = 'Ms' # what the?
    multiplier_values = np.asarray([[1.0]])
    lat = np.asarray([[-24.0]])
    lon = np.asarray([[144.0]])
    nc_name = 'example.nc'
    save_multiplier(multiplier_name, multiplier_values, lat, lon, nc_name)

# Get rid of this.  It's moving to processMultipliers
def course_yasi_img():
    tl_y = np.asarray([-16]) # top left y
    tl_x = np.asarray([140]) # top left x
    dx = 2
    dy = -2
    multiplier_values = np.zeros((2, 3))

    for index in indices:
        multiplier_values.fill(index)
        img_name = 'm4_' + indices[index]['dir'] + '.img'
        pM.createRaster(multiplier_values, tl_x, tl_y,
                                 dx, dy,
                                 filename=img_name)

def course_yasi_nc():
    """
    This was the wrong file format.
    :return:
    """
    multiplier_name = 'Ms' # what the?
    lat = np.asarray([ -23, -20, -17, -14, -11, -8, -5])
    lon = np.asarray([137, 140, 143, 146, 149, 152, 155, 158])
    multiplier_values = np.zeros(([lat.shape[0], lon.shape[0]]))

    for index in indices:
        multiplier_values.fill(index)
        nc_name = 'syn_' + indices[index]['dir'] + '.nc'
        save_multiplier(multiplier_name, multiplier_values, lat, lon, nc_name)

# -------------------------------------------------------------
if __name__ == "__main__":
    #example()
    #course_yasi_nc()
    course_yasi_img()