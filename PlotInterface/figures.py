import sys
import numpy as np
import scipy.stats as stats
import wind.windmodels as windmodels
from matplotlib.figure import Figure


class WindProfileFigure(Figure):

    def __init__(self, lat, lon, eP, cP, rMax, beta, beta1=1.5, beta2=1.4):
        Figure.__init__(self)
        self.R = np.array(range(1, 201), 'f')
        self.lat = lat
        self.lon = lon
        self.rMax = rMax
        self.eP = eP
        self.cP = cP
        self.beta = beta
        self.beta1 = beta1
        self.beta2 = beta2

    def plot(self, profileType=None):
        profiles = []

        if profileType:
            profiles.append(profileType)
        else:
            for p in windmodels.PROFILES:
                profiles.append(p)

        ax = self.add_subplot(1, 1, 1)
        ax.hold(True)
        legend = []

        for name in profiles:
            try:
                cls = windmodels.profile(name)
                params = windmodels.profileParams(name)
                values = [getattr(self, p) for p in params if hasattr(self, p)]
                profile = cls(self.lat, self.lon, self.eP, self.cP,
                              self.rMax, *values)
                V = profile.velocity(self.R)
                ax.plot(self.R, abs(V), linewidth=2)
                legend.append(name.capitalize())
            except TypeError:
                pass

        ax.legend(legend)
        ax.grid()
        ax.set_xlabel('Radius (km)', fontsize=14)
        ax.set_ylabel('Wind speed (m/s)', fontsize=14)
        ax.set_title((r'$P_c = %d\hspace{0.5}hPa,\hspace{1} P_e' +
                      r'= %d \hspace{0.5} hPa,\hspace{1} R_{max}' +
                      r'= %d \hspace{0.5}km$') %
                    (self.cP / 100., self.eP / 100., self.rMax))


class RegressionFigure(Figure):

    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []

    def add(self, data, label='x', title='x', transform=lambda x: x):
        self.subfigures.append((data, label, title, transform))

    def prepareData(self, data, transform):
        z = data
        i = np.where((z[0, :] < sys.maxint) & (z[1, :] < sys.maxint))[0]
        z = transform(z[:, i])
        return z

    def subplot(self, axes, data, label='x', title='x', transform=lambda x: x):
        color = axes._get_lines.color_cycle

        z = self.prepareData(data, transform)

        k = color.next()
        axes.plot(z[0,:], z[1,:], '.', alpha=0.8, color=k)

        def filterOutliers(data, m=12.):
            d = np.abs(data - np.median(data))
            mdev = np.median(d)
            s  = d/mdev if mdev else 0
            return data[s<m]

        xmin = filterOutliers(z[0,:]).min()
        xmax = filterOutliers(z[0,:]).max()

        m, c, r, p, e = stats.linregress(z)
        x = np.array([xmin, xmax])
        y = m * x + c

        k = color.next()
        axes.plot(x, x, '-', color=k)

        k = color.next()
        axes.plot(x, y, '-', color=k, label='r = %5.3f' % r)

        axes.set_xlim(xmin, xmax)
        axes.set_ylim(xmin, xmax)
        axes.set_xticks(np.linspace(xmin, xmax, 7))
        axes.set_yticks(np.linspace(xmin, xmax, 7))

        axes.set_ylabel('$' + label + '(t)$')
        axes.set_xlabel('$' + label + '(t-1)$')
        axes.set_title(title)

        axes.set_aspect('equal')

        legend = axes.legend(loc=2)
        legend.get_frame().set_alpha(0.5)

    def plot(self):
        n = len(self.subfigures)
        w,h = self.get_size_inches()
        self.set_size_inches(w, n*h)
        for i, subfigure in enumerate(reversed(self.subfigures)):
            axes = self.add_subplot(n, 1, i)
            self.subplot(axes, *subfigure)


class LaggedRegressionFigure(RegressionFigure):

    def prepareData(self, data, transform):
        z = np.array([data[1:], data[:-1]])
        i = np.where((z[0, :] < sys.maxint) & (z[1, :] < sys.maxint))[0]
        z = transform(z[:, i])
        return z


class PressureFigure(LaggedRegressionFigure):

    def plot(self, pressures, pressureRates):
        self.add(pressures, 'p', 'Pressure')
        self.add(pressureRates, r'\Delta p', 'Pressure rate of change')
        super(PressureFigure, self).plot()


class SpeedFigure(LaggedRegressionFigure):

    def plot(self, speeds, speedRates):
        self.add(speeds, 'v', 'Speed')
        self.add(speedRates, r'\Delta v', 'Speed rate of change (Acceleration)')
        super(SpeedFigure, self).plot()


class BearingFigure(LaggedRegressionFigure):

    def plot(self, bearings, bearingRates):
        def transform(z):
            return np.cos(np.radians(z))
        self.add(bearings, r'\cos(\theta)', 'Bearing', transform)
        self.add(bearingRates, r'\Delta\theta',
                 'Bearing rate of change')
        super(BearingFigure, self).plot()


class FrequencyFigure(RegressionFigure):

    def plot(self, years, frequency):
        data = np.array([years, frequency])
        self.add(data, '', 'Annual Frequency')
        super(FrequencyFigure, self).plot()


def saveFigure(figure, filename):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    canvas = FigureCanvas(figure)
    canvas.print_figure(filename)


def saveWindProfilesFigure(lat, lon, eP, cP, rMax, beta,
                           filename='docs/windprofiles.png'):
    fig = WindProfileFigure(lat, lon, eP, cP, rMax, beta)
    fig.plot()
    saveFigure(fig, filename)


def savePressureFigure(pressures, pressureRates, filename='docs/prs_corr.png'):
    fig = PressureFigure()
    fig.plot(pressures, pressureRates)
    saveFigure(fig, filename)


def saveSpeedFigure(pressures, pressureRates, filename='docs/spd_corr.png'):
    fig = SpeedFigure()
    fig.plot(pressures, pressureRates)
    saveFigure(fig, filename)


def saveBearingFigure(bearings, bearingRates, filename='docs/bear_corr.png'):
    fig = BearingFigure()
    fig.plot(bearings, bearingRates)
    saveFigure(fig, filename)


def saveFrequencyFigure(years, frequency, filename='docs/bear_corr.png'):
    fig = FrequencyFigure()
    fig.plot(years, frequency)
    saveFigure(fig, filename)


def main():
    import Utilities.files as files
    from os.path import join as pjoin, dirname, normpath
    baseDir = normpath(pjoin(dirname(__file__), '..'))
    inputPath = pjoin(baseDir, 'output', 'process')

    saveWindProfilesFigure(-12., 130., 100700., 95000., 30., 1.6)

    pressureRates = files.flLoadFile(pjoin(inputPath, 'pressure_rate'))
    pressures = files.flLoadFile(pjoin(inputPath, 'all_pressure'))
    savePressureFigure(pressures, pressureRates)

    speedRates = files.flLoadFile(pjoin(inputPath, 'speed_rate'))
    speeds = files.flLoadFile(pjoin(inputPath, 'all_speed'))
    saveSpeedFigure(speeds, speedRates)

    bearingRates = files.flLoadFile(pjoin(inputPath, 'bearing_rate'))
    bearings = files.flLoadFile(pjoin(inputPath, 'all_bearing'))
    saveBearingFigure(bearings, bearingRates)

if __name__ == "__main__":
    main()
