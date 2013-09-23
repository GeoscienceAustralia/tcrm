import numpy as np
from Utilities.config import ConfigParser

def colReadCSV(configFile, dataFile, source):
    """
    Loads a csv file containing 'column' data into a record
    (numpy) array with columns labelled by 'fields'.
    """
    config = ConfigParser()
    config.read(configFile)
    delimiter = config.get(source, 'FieldDelimiter')
    numHeadingLines = config.getint(source, 'NumberOfHeadingLines')
    cols = config.get(source, 'Columns').split(delimiter)

    usecols = [i for i,c in enumerate(cols) if c != 'skip']

    data = np.genfromtxt(dataFile, dtype=None, delimiter=delimiter,
            usecols=usecols, comments=None, skip_header=numHeadingLines, 
            autostrip=True)

    data.dtype.names = [c for c in cols if c != 'skip']

    return data
