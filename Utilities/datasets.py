from urllib2 import urlopen
from cStringIO import StringIO
from Utilities.config import ConfigParser
from os.path import isfile, splitext, join as pjoin

DATASETS = []


class DataSet(object):

    def __init__(self, name, url, path):
        self.name = name
        self.url = url
        self.path = path
        self.filename = pjoin(path, url.split('/')[-1])
        self.compression = None

        base, ext = splitext(self.filename)
        if ext in ['.gz', '.zip']:
            self.compression = ext[1:]
            self.filename = base

    def download(self, callback=None):
        if self.isDownloaded():
            return

        urlfile = urlopen(self.url, timeout=3)
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

        with open(self.filename, 'wb') as outfile:
            outfile.write(data.read())

    def isDownloaded(self):
        return isfile(self.filename)



def loadDatasets():
    config = ConfigParser()
    datasets = config.get('Input', 'Datasets').split(',')
    global DATASETS
    for dataset in datasets:
        url = config.get(dataset, 'URL')
        path = config.get(dataset, 'location')
        data = DataSet(dataset, url, path)
        DATASETS.append(data)


def checkAndDownload(callback=None):
    for dataset in DATASETS:
        dataset.download()


loadDatasets()


if __name__ == '__main__':
    def status(fn, done, size):
        status = r"%s  %10d  [%3.2f%%]" % (fn, done, done * 100. / size)
        status = status + chr(8)*(len(status)+1)
        print status,
    checkAndDownload(status)
