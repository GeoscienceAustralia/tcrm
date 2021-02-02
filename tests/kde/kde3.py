import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

import Utilities.KPDF as KPDF
import statsmodels.nonparametric.kernel_density as smkde
import statsmodels.nonparametric.bandwidths as smbw

import seaborn as sns
sns.set_context("paper")

xrvs = [stats.norm(loc=1.0, scale=0.2), stats.norm(loc=3.0, scale=0.1)]
yrvs = [stats.norm(loc=0.0, scale=0.2), stats.norm(loc=1.0, scale=0.3)]

N = 3000
M = N / (len(xrvs) + len(yrvs))

rvs = []
for xrv in xrvs:
    for yrv in yrvs:
        rvs.append(np.hstack([
            xrv.rvs(size=(M, 1)),
            yrv.rvs(size=(M, 1))]))
rvs = np.vstack(rvs)


def pdf(x, y):
    return (sum([rv.pdf(x) for rv in xrvs]) *
            sum([rv.pdf(y) for rv in yrvs]))


def plot(ax, Z, title):
    ax.scatter(rvs[:, 0], rvs[:, 1], alpha=0.3, s=1, color='black')
    ax.imshow(Z, aspect=1, origin='lower',
              cmap=sns.light_palette((210, 90, 60),
                                     input='husl',
                                     as_cmap=True),
              extent=(rvs[:, 0].min(),
                      rvs[:, 0].max(),
                      rvs[:, 1].min(),
                      rvs[:, 1].max()))
    ax.contour(x, y, Z, linewidth=1)
    ax.set_title(title)

x_flat = np.r_[rvs[:, 0].min():rvs[:, 0].max():128j] #pylint: disable=E1127
y_flat = np.r_[rvs[:, 1].min():rvs[:, 1].max():128j] #pylint: disable=E1127
x, y = np.meshgrid(x_flat, y_flat)
grid = np.hstack([x.reshape(-1, 1), y.reshape(-1, 1)])

n, d = rvs.shape
bw_kpdf = KPDF.MPDFOptimumBandwidth(rvs)
bw_scott = n ** (-1.0 / (d + 4))
bw_silverman = (n * (d + 2) / 4.) ** (-1. / (d + 4))

scale = np.array([[x_flat.ptp(), y_flat.ptp()]])
print(('scale: %s' % scale))

print(('bandwidth kpdf=%f scott=%f silverman=%f' %
      (bw_kpdf, bw_scott, bw_silverman)))

fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(8, 11))
axes = axes.flat

p = pdf(x, y)
p = p / p.max()
plot(axes[0], p, 'True density')

kde = stats.kde.gaussian_kde(rvs.T)
print((kde.covariance))
z = kde(grid.T)
z = z.reshape(x.shape) / z.max()
plot(axes[1], z, 'Scipy ($\ell_2$ norm: %.3f)' % np.linalg.norm((p - z).flat))

bw = bw_kpdf
# bw = bw_scott

w = KPDF.MPDFGaussian(rvs, grid, bw_kpdf)
w = w.reshape(x.shape) / w.max()
plot(axes[2], w, 'KPDF bw:kpdf ($\ell_2$ norm: %.3f)' % np.linalg.norm((p -
     w).flat))

w = KPDF.MPDFGaussian(rvs, grid, bw_kpdf / 2)
w = w.reshape(x.shape) / w.max()
plot(axes[3], w, 'KPDF bw:kpdf/2 ($\ell_2$ norm: %.3f)' % np.linalg.norm((
    p - w).flat))

w = KPDF.MPDFGaussian(rvs, grid, bw_scott)
w = w.reshape(x.shape) / w.max()
plot(axes[4], w, 'KPDF bw:scott ($\ell_2$ norm: %.3f)' % np.linalg.norm((p -
     w).flat))

w = KPDF.MPDFGaussian(rvs, grid, bw_scott / 2)
w = w.reshape(x.shape) / w.max()
plot(axes[5], w, 'KPDF bw:scott/2 ($\ell_2$ norm: %.3f)' % np.linalg.norm((
    p - w).flat))

dens = smkde.KDEMultivariate(rvs, 'cc', bw='cv_ml')
print("SM bandwidth (cv_ml): " + repr(dens.bw))
w = dens.pdf(grid)
w = w.reshape(x.shape) / w.max()
plot(axes[6], w, 'SM bw:CVML ($\ell_2$ norm: %.3f)' % np.linalg.norm((p -
     w).flat))

dens = smkde.KDEMultivariate(rvs, 'cc', bw='cv_ls')
print("SM bandwidth (cv_ls): " + repr(dens.bw))
w = dens.pdf(grid)
w = w.reshape(x.shape) / w.max()
plot(axes[7], w, 'SM bw:CVLS ($\ell_2$ norm: %.3f)' % np.linalg.norm((p -
     w).flat))

plt.savefig('fig3.png')
plt.show()
