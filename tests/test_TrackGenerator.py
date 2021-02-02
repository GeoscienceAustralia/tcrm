import unittest
from numpy.testing import *

# from TrackGenerator import TrackGenerator


class TestTrackGenerator(unittest.TestCase):

    def testGeneratePath(self):
        assert_almost_equal(list(range(10)), list(range(10)))
        pass

if __name__ == "__main__":
    suite = unittest.makeSuite(TestTrackGenerator, 'test')
    unittest.TextTestRunner().run(suite)
