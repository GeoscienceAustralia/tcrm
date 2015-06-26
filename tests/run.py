import unittest
import logging

if 'discover' not in dir(unittest.TestLoader):
    import py26compat
    unittest.TestLoader = py26compat.TestLoader

logging.disable(logging.CRITICAL)

if __name__ == '__main__':
    suite = unittest.TestLoader().discover('.')
    unittest.TextTestRunner(verbosity=2).run(suite)
