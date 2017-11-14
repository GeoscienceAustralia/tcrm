__author__ = 'duncan'

from os.path import join as pjoin, exists
from ProcessMultipliers import processMultipliers as pM


def gen_one_syns():
    tl_y = -20   # top left y
    tl_x = 138   # top left x
    delta = 2
    dir_path = '.'
    shape=(2, 4)

    pM.generate_syn_mult_img(tl_x, tl_y, delta, dir_path, shape,
                          every_fill=1.0)

if __name__ == "__main__":
    gen_one_syns()
