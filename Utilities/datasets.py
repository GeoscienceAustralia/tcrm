from urllib2 import urlopen, URLError
from cStringIO import StringIO
from Utilities.config import ConfigParser
from os.path import isfile, splitext, join as pjoin

DATASETS = []


class DataSet(object):

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
        if self.isDownloaded():
            return

        try:
            urlfile = urlopen(self.url, timeout=5)
            meta = urlfile.info()
            data = StringIO()

            size = int(meta.getheaders('Content-Length')[0])
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
        return isfile(pjoin(self.path, self.filename))


def loadDatasets():
    config = ConfigParser()
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
    for dataset in DATASETS:
        dataset.download(callback)


loadDatasets()


if __name__ == '__main__':
    def status(fn, done, size):
        status = r"%s  %10d  [%3.2f%%]" % (fn, done, done * 100. / size)
        status = status + chr(8)*(len(status)+1)
        print status,
    checkAndDownload(status)
