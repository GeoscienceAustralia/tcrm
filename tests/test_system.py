import unittest
from unittest.mock import patch
import os.path
import tempfile
import imageio
import numpy as np

import Utilities.config
import Utilities.pathLocator
import Evaluate.interpolateTracks
import tcevent

def decimate(factor):
    """
    Patch to reduce number of time steps in simulation.

    Ideally instead there should be a config option to change the time step,
    and that supports time steps larger than the observation intervals.

    Instead, this works by downsampling as soon as a track is interpolated,
    ie. thinning out the timeseries that will be input for simulation.

    Incidentally, for large factors the output plot contains the aesthetic
    effect of superimposing a sequence of snapshots of the cyclone, revealing
    wind model internal structure that is not apparent in the full aggregate.
    """
    raw_interpolate = Evaluate.interpolateTracks.interpolate
    def wrapper(*args, **kwargs):
        result = raw_interpolate(*args, **kwargs)
        result.data = result.data[::factor]
        return result
    return patch('Evaluate.interpolateTracks.interpolate', wrapper)

class Yasi_Example(unittest.TestCase):
    def setUp(self):
        # Prevent interference with other tests
        Utilities.config.reset()
        self.addCleanup(Utilities.config.reset)

        self.tmpdir = tempfile.TemporaryDirectory()
        #self.addCleanup(self.tmpdir.cleanup)

        self.configFile = os.path.join(
                            Utilities.pathLocator.getRootDirectory(),
                            "example/yasi.ini")

        # Use singleton property to propagate config updates to other modules
        config = Utilities.config.ConfigParser()
        config.read(self.configFile)
        config['Output']['Path'] = self.tmpdir.name
        config['WindfieldInterface']['PlotOutput'] = 'True'

    @decimate(100)
    def test_scenario(self):
        fname = os.path.join(self.tmpdir.name, "plots/maxwind.png")

        tcevent.main(self.configFile) # run simulation

        with self.subTest("Reading plots/maxwind.png"):
            img = imageio.imread(fname)

        with self.subTest("Inspecting plot (after track decimation)"):
            x, y, bands = img.shape
            pixels = x * y
            self.assertGreater(x, 50) # check image has extent
            self.assertGreater(y, 50)
            black = np.all(img[:,:,:3] == [0] * 3, axis=-1)
            white = np.all(img[:,:,:3] == [255] * 3, axis=-1)
            color = ~(black | white)
            self.assertGreater(black.sum(), 100) # linework exists
            self.assertLess(black.sum() / pixels, 0.2) # but doesn't dominate
            self.assertGreater(white.sum() / pixels, 0.2) # substantial space
            self.assertGreater(color.sum() / pixels, 0.05) # significant color

        from time import sleep
        sleep(1)
        
if __name__ == '__main__':
    unittest.main()
