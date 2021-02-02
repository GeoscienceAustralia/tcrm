import numpy as np
import wind.windmodels as windmodels
from Utilities.maputils import ModelGrid
from matplotlib.figure import Figure
import seaborn as sns
sns.set(style="ticks")
pal = sns.color_palette('Paired', 10)
sns.set_palette(pal)
import pdb

class WindProfileFigure(Figure):

    def __init__(self, lat, lon, eP, cP, rMax, beta, dcP=200., vFm=2.5, 
                 rGale=250000., beta1=1.5, beta2=1.4):
        Figure.__init__(self)
        self.R = np.array(list(range(1, 201)), 'f')
        self.lat = lat
        self.lon = lon
        self.rMax = rMax
        self.eP = eP
        self.cP = cP
        self.beta = beta
        self.dcP = dcP
        self.vFm = vFm
        self.rGale = rGale
        self.beta1 = beta1
        self.beta2 = beta2
        self.grid = ModelGrid(self.lon, self.lat, 2.0, 0.01)
        self.grid.makeGrid()

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
        nx = self.grid.R.shape[0] 
        idx = (nx - 1)/2
        for name in profiles:
            try:
                print name
                cls = windmodels.profile(name)
                params = windmodels.profileParams(name)
                values = [getattr(self, p) for p in params if hasattr(self, p)]
                profile = cls(self.grid, self.eP, self.cP,
                              self.rMax, *values)
                V = profile.velocity(self.grid)
                ax.plot(self.grid.R[idx,:]/1000., abs(V[idx,:]), linewidth=2)
                legend.append(name.capitalize())
            except TypeError:
                pass
            except:
                raise

        ax.legend(legend)
        ax.set_xlabel('Radius (km)', fontsize=14)
        ax.set_ylabel('Wind speed (m/s)', fontsize=14)
        ax.set_title((r'$P_c = %d\hspace{0.5}hPa,\hspace{1} P_e' +
                      r'= %d \hspace{0.5} hPa,\hspace{1} R_{max}' +
                      r'= %d \hspace{0.5}km$') %
                      (self.cP/100., self.eP/100., self.rMax/1000.))
        ax.set_xlim((0, 200))

def main():
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from os.path import join as pjoin, dirname, normpath

    lat = -12.
    lon = 130.
    rMax = 20000.
    eP = 100700.
    cP = 95000.
    beta = 1.6
    dcP = 200.
    vFm = 2.5
    rGale = 100000.

    filename = 'windprofiles.png'
    path = pjoin(normpath(pjoin(dirname(__file__), '..', 'docs')), filename)

    fig = WindProfileFigure(lat, lon, eP, cP, rMax, beta, dcP, vFm, rGale)
    canvas = FigureCanvas(fig)

    fig.plot()

    print(('saving wind profiles figure to %s' % path))
    canvas.print_figure(path)

if __name__ == "__main__":
    main()
