import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from mpl_toolkits.basemap import Basemap
import Utilities.KPDF as KPDF

ff = os.path.join('..', 'test_data', 'kde_origin_lonLat.pkl')
lonlat = pickle.load(open(ff, 'rb'))

xmin = 70.0
xmax = 180.0
ymin = -36.
ymax = 0.

fig, axes = plt.subplots(nrows=2, ncols=1)

for ax in axes.flat:
    m = Basemap(ax=ax, projection='cyl', resolution='i', llcrnrlon=xmin,
                urcrnrlon=xmax, llcrnrlat=ymin, urcrnrlat=ymax)
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(fill_color='#ffffff')
    m.fillcontinents(color='#dedcd2', lake_color='#ffffff')
    meridians = np.arange(xmin, xmax, 10.)
    parallels = np.arange(ymin, ymax, 10.)
    m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=9, linewidth=0.2)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=9, linewidth=0.2)
    ax.grid(True)

xy = lonlat
n, d = xy.shape
bw_kpdf = KPDF.MPDFOptimumBandwidth(xy)
bw_scott = n ** (-1.0 / (d + 4))
bw_silverman = (n * (d + 2) / 4.) ** (-1. / (d + 4))

X, Y = np.mgrid[xmin:xmax:30j, ymin:ymax:30j]
positions = np.vstack([X.ravel(), Y.ravel()]).T

# KDE

z = KPDF.MPDFGaussian(xy, positions, bw_kpdf)
z = np.reshape(z, X.shape)
print(('z max = %f' % z.max()))

ax = axes[0]
ax.imshow(np.rot90(z) / z.max(), cmap=plt.cm.PuRd,
          extent=[xmin, xmax, ymin, ymax])
ax.scatter(lonlat[:, 0], lonlat[:, 1], s=2, marker='o',
           facecolor='darkred', edgecolor='darkred', alpha=0.8)
ax.contour(X, Y, np.rot90(z) / z.max(), colour='magenta')
ax.set_title('KPDF')

# Scipy

kde = stats.gaussian_kde(lonlat.T)
Z = np.reshape(kde(positions.T), X.shape)
print(('Z max = %f' % Z.max()))

ax = axes[1]
ax.imshow(np.rot90(Z) / Z.max(), cmap=plt.cm.PuRd,
          extent=[xmin, xmax, ymin, ymax])
ax.scatter(lonlat[:, 0], lonlat[:, 1], s=2, marker='o',
           facecolor='darkred', edgecolor='darkred', alpha=0.8)
ax.contour(X, Y, np.rot90(Z) / Z.max(), colour='magenta')
ax.set_title('Scipy')

plt.savefig('fig2.pdf')
plt.show()
