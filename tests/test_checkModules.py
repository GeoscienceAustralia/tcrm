import unittest

from tcrm import checkModules

class TestCheckModules(unittest.TestCase):

    def test_missingModule(self):
        """checkModules returns False for missing module"""
        self.assertFalse(checkModules(['cmaps']))

    def test_knownModule(self):
        """checkModules returns true if all modules available"""
        self.assertTrue(checkModules(['os']))

    def test_notListInput(self):
        """checkModules raises TypeError if not list input"""
        self.assertRaises(TypeError, checkModules, 'os')

if __name__ == "__main__":
    unittest.main()
