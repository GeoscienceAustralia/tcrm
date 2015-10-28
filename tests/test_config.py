import unittest

import io

import config

class TestParsers(unittest.TestCase):

    def setUp(self):
        self.inputLine = "A,B,C,D"
        self.outputList = ["A", "B", "C", "D"]
        self.singleInputLine = "A"
        self.inputList = ['X', 'Y', 'Z']
        self.outputLine = 'X,Y,Z'

    def test_parseBool(self):
        """Parse boolean text strings"""
        self.assertTrue(config.parseBool("True"))
        self.assertFalse(config.parseBool("true"))
        self.assertFalse(config.parseBool("10"))
        self.assertFalse(config.parseBool("False"))

    def test_parseList(self):
        """Parse comma-separated line into a list"""
        self.assertEquals(config.parseList(self.inputLine), self.outputList)
        self.assertEquals(config.parseList(self.singleInputLine), ["A"])

    def test_formatList(self):
        """Test conversion to comma-joined string"""
        self.assertEquals(config.formatList(self.inputList), self.outputLine)

class TestConfigParser(unittest.TestCase):

    def setUp(self):
        self.testConfig = """
[Boolean]
Option1=True
Option2=False

[Integer]
Option1=1
Option2=10

[Float]
Option1=0.5
Option2=2.5

[String]
Option=randomstring

[List]
Option1=A,B,C,D,E
Option2=A

[Evaluate]
Option={'x':1.0,'y':0.5}
"""
        self.config = config.ConfigParser(defaults=self.testConfig)
        self.evalopt = {'x':1.0,'y':0.5}
        self.listopt1 = ["A", "B", "C", "D", "E"]
        self.listopt2 = ["A"]

    def test_geteval(self):
        """Test geteval returns correct object"""
        self.assertIsInstance(self.config.geteval('Evaluate', 'Option'), dict)
        self.assertEqual(self.config.geteval('Evaluate', 'Option'), self.evalopt)

#    def test_items(self):
#        """Test ConfigParser.items method"""
#        self.assertEqual(self.config.items('Float'), [('option1', 0.5), ('option2', 2.5)])

class TestOldStyleConfig(unittest.TestCase):
    def setUp(self):
        self.testConfig = """
[Boolean]
Option1=True
Option2=False

[Integer]
Option1=1
Option2=10

[Float]
Option1=0.5
Option2=2.5

[String]
Option=randomstring

[List]
Option1=A,B,C,D,E
Option2=A

[Evaluate]
Option={'x':1.0,'y':0.5}
"""
        self.configFile = io.BytesIO(self.testConfig)

    def test_cnfReadIniValue(self):
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Integer', 'Option1'), 1)
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Integer', 'Option2'), 10)
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Boolean', 'Option1'), True)
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Boolean', 'Option2'), False)
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Float', 'Option1'), 0.5)
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Float', 'Option2'), 2.5)
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'List', 'Option1'), "A,B,C,D,E")
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'List', 'Option2'), "A")
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'String', 'Option'), 'randomstring')
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Evaluate', 'Option'), {'x':1.0,'y':0.5})

        
    def test_returnDefault(self):
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'Integer', 'Option3', 1), 1)
        self.assertEqual(config.cnfGetIniValue(self.configFile, 'String', 'Option2', 'otherstring'), 'otherstring')

if __name__ == "__main__":
    unittest.main()
