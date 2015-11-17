"""
:mod:`dynarray` -- a dynamic record array
================================================

.. module:: dynarray
    :synopsis: A dynamic record array to enable appending
               and extending an existing record array with
               additional records.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import numpy as np

class DynamicRecArray(object):
    """
    A dynamic record array to enable appending/extending an array
    """
    def __init__(self, dtype):
        self.dtype = np.dtype(dtype)
        self.length = 0
        self.size = 10
        self._data = np.empty(self.size, self.dtype)

    def __len__(self):
        return self.length

    def append(self, rec):
        """
        Append a record to the array.

        :param rec: new record to append to the existing array

        """
        if self.length == self.size:
            self.size = int(1.5 * self.size)
            self._data = np.resize(self._data, self.size)
        self._data[self.length] = rec
        self.length += 1

    def extend(self, recs):
        """
        Extend a record array  with many records.

        :param recs: list of new records to append to the array.

        """
        for rec in recs:
            self.append(rec)

    @property
    def data(self):
        """
        Convenience function to return the data stored in the array.
        """
        return self._data[:self.length]

