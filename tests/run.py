import unittest
import logging

logging.disable(logging.CRITICAL)

if __name__ == '__main__':
    suite = unittest.TestLoader().discover('.')
    unittest.TextTestRunner(verbosity=2).run(suite)
