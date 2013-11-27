from urllib2 import urlopen
from cStringIO import StringIO

class DataSet(object):

    def __init__(self, name, url):
        self.name = name
        self.url = url

    def download(self, callback=None):
        fn = self.url.split('/')[-1]
        raw = StringIO()
        u = urlopen(self.url)
        meta = u.info()
        if callback is None:
            def callback(x, y):
                pass
        size = int(meta.getheaders('Content-Length')[0])
        done = 0
        while True:
            buf = u.read(8192)
            if not buf:
                break
            raw.write(buf)
            done += len(buf)
            callback(done, size)
        if fn.endswith('gz'):
            from gzip import GzipFile as zipfile
            raw.seek(0)
            with open(fn[:-3], 'wb') as outfile:
                outfile.write(zipfile(fileobj=raw).read())

IBTRACS = DataSet('IBTRACS', 'ftp://eclipse.ncdc.noaa.gov/pub/ibtracs/v03r05/wmo/csv/Allstorms.ibtracs_wmo.v03r05.csv.gz')

def status(done, size):
    status = r"%10d  [%3.2f%%]" % (done, done * 100. / size)
    status = status + chr(8)*(len(status)+1)
    print status,

IBTRACS.download(status)
