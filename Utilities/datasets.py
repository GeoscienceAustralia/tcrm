"""
:mod:`datasets` -- manage downloaded datasets
=============================================

"""
from urllib.request import urlopen
from urllib.error import URLError
from io import BytesIO
from Utilities.config import ConfigParser
from os.path import isfile, splitext, join as pjoin

DATASETS = []


class DataSet(object):
    """
    Download a dataset from a url, store in a specified path with a
    given filename.

    :param str name: Name of the dataset to be downloaded.
    :param str url: URL of the dataset.
    :param str path: path name for the storage location of the dataset
    :param filename: name of the file to be saved (can be different from
                     the name of the dataset).
    :type filename: str or None

    """
    def __init__(self, name, url, path, filename=None):
        self.name = name
        self.url = url
        self.path = path
        self.compression = None

        base, ext = splitext(url.split('/')[-1])
        if ext in ['.gz', '.zip']:
            self.compression = ext[1:]
            self.filename = filename or base
        else:
            self.filename = filename or url.split('/')[-1]

    def download(self, callback=None):
        """
        Execute the download of the dataset. If the dataset
        has already been downloaded (the :attr:`isDownloaded` will be set
        to True), then don't try to download.

        :param callback: Callback function (for reporting status to STDOUT).
        :type callback: function

        :raises IOError: Unable to download the dataset.

        """

        if self.isDownloaded():
            return

        try:
            urlfile = urlopen(self.url, timeout=5)
            meta = urlfile.info()
            data = BytesIO()

            size = int(meta['Content-Length'])
            done = 0
            while True:
                buf = urlfile.read(8192)
                if not buf:
                    break
                data.write(buf)
                done += len(buf)
                if callback is not None:
                    callback(self.filename, done, size)

            data.seek(0)

            if self.compression == 'gz':
                from gzip import GzipFile
                data = GzipFile(fileobj=data)

            if self.compression == 'zip':
                from zipfile import ZipFile
                data = ZipFile(data)

            with open(pjoin(self.path, self.filename), 'wb') as outfile:
                outfile.write(data.read())

        except URLError:
            raise IOError('Cannot download file')

    def isDownloaded(self):
        """
        Determine if a file has already been downloaded

        :returns: `True` if the file exists, `False` otherwise.

        """
        return isfile(pjoin(self.path, self.filename))


def loadDatasets(configFile):
    """
    Load the details of the datasets to be downloaded from the
    configuration settings. This updates the :data:`DATASETS`
    list.

    """

    config = ConfigParser()
    config.read(configFile)
    datasets = config.get('Input', 'Datasets').split(',')

    global DATASETS
    for dataset in datasets:
        url = config.get(dataset, 'URL')
        path = config.get(dataset, 'path')
        if config.has_option(dataset, 'filename'):
            filename = config.get(dataset, 'filename')
        else:
            filename = None
        data = DataSet(dataset, url, path, filename)
        DATASETS.append(data)


def checkAndDownload(callback=None):
    """
    Check the :data:`DATASETS` list and download each, if it
    has not previously been downoladed.

    :param callback: Callback function (for reporting status to STDOUT).
    :type callback: function


    """

    for dataset in DATASETS:
        dataset.download(callback)


if __name__ == '__main__':
    def status(fn, done, size):
        status = r"%s  %10d  [%3.2f%%]" % (fn, done, done * 100. / size)
        status = status + chr(8)*(len(status)+1)
        print(status, end=' ')
    checkAndDownload(status)
