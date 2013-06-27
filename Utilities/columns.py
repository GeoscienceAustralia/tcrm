import logging
import numpy as np
from config import ConfigParser

log = logging.getLogger()

def colReadCSV(configFile, dataFile, source):
    """
    Loads a csv file containing 'column' data into a dictionary of
    (numpy) arrays with keys labelled by 'fields'. Values of the
    dictionary are arrays containing the data in the columns in the
    file.  This is necessary as the csv.reader() function
    (and hence DictReader) returns all values as strings. No automatic
    data type conversion is performed.
    Input: dataFile - path to the csv file to read
           source - name of a data source format describing the data.
                    The source should have a corresponding section in
                    gConfigFile.
    Output: data - a dictionary of arrays, with keys given by the
                   column names specified in the configuration file.
    Example: data = colReadCSV(dataFile, source)
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
