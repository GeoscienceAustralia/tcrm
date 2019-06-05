import matplotlib as mpl
mpl.rcParams['legend.fancybox'] = True

from matplotlib.figure import Figure
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, date2num
from matplotlib.ticker import MultipleLocator



class TimeSeriesFigure(Figure):
    """
    Base class for plotting time series data
    """
    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []
        self.subplots_adjust(bottom=0.1, right=0.8)
        self.dayLocator = DayLocator()
        self.hourLocator = HourLocator(interval=3)
        self.dateFormat = DateFormatter('%H:%M \n %Y-%m-%d')
        self.set_size_inches(12, 4)
        self.dirTickLocator = MultipleLocator(45)


    def make_patch_spines_invisible(self, axes):
        axes.set_frame_on(True)
        axes.patch.set_visible(False)
        for sp in axes.spines.values():
            sp.set_visible(False)

    def add(self, dt, ydata, yrange, ylabel, title):
        self.subfigures.append((dt, ydata, yrange, ylabel, title))

    def subplot(self, axes, subfigure, color):

        dt, ydata, yrange, ylabel, title = subfigure
        dt = date2num(dt)

        ymin, ymax = yrange

        axes.set_ylim(ymin, ymax)
        p = axes.plot(dt, ydata, '-', color=color, label=ylabel, linewidth=2)

        self.formatXticks(axes)
        self.addGrid(axes)
        axes.set_ylabel(ylabel)
        axes.yaxis.label.set_color(p[0].get_color())
        axes.tick_params(axis='y', colors=p[0].get_color())
        axes.tick_params(axis='x', colors=p[0].get_color())
        return p[0]

    def subplot_twinx(self, axes, subfigure, color, position):

        axes.spines['right'].set_position(("axes", position))
        self.make_patch_spines_invisible(axes)
        axes.spines['right'].set_visible(True)

        dt, ydata, yrange, ylabel, title = subfigure
        dt = date2num(dt)

        ymin, ymax = yrange

        p = axes.plot(dt, ydata, '-', color=color, label=ylabel, linewidth=2)
        axes.set_ylabel(ylabel)
        axes.set_ylim(ymin, ymax)
        axes.yaxis.label.set_color(p[0].get_color())
        axes.tick_params(axis='y', colors=p[0].get_color())
        return p[0]

    def plot(self):

        axes = self.add_subplot(111)
        pc = mpl.rcParams['axes.prop_cycle']
        color = iter(pc.by_key()['color'])
        ax1 = axes.twinx()
        ax2 = axes.twinx()


        ax2.spines['right'].set_position(("axes", 1.1))
        self.make_patch_spines_invisible(ax2)
        ax2.spines['right'].set_visible(True)


        dt, ydata0, yrange0, ylabel0, title0 = self.subfigures[0]
        dt, ydata1, yrange1, ylabel1, title1 = self.subfigures[1]
        dt, ydata2, yrange2, ylabel2, title2 = self.subfigures[2]
        dt = date2num(dt)

        p0, = axes.plot(dt, ydata0, color=next(color),
                        label=title0, linewidth=2, alpha=0.5)
        p1, = ax1.plot(dt, ydata1, color=next(color),
                        label=title1, linewidth=2, alpha=0.5)
        p2, = ax2.plot(dt, ydata2, color=next(color),
                        label=title2, linewidth=2, alpha=0.5)

        axes.set_xlabel('Time')
        axes.set_ylabel(ylabel0)
        axes.set_ylim(900, 1020)
        ax1.set_ylabel(ylabel1)
        ax1.set_ylim(0, 100)
        ax2.set_ylabel(ylabel2)
        ax2.set_ylim(0, 360)

        axes.yaxis.label.set_color(p0.get_color())
        ax1.yaxis.label.set_color(p1.get_color())
        ax2.yaxis.label.set_color(p2.get_color())
        ax2.yaxis.set_major_locator(self.dirTickLocator)

        axes.tick_params(axis='y', colors=p0.get_color())
        axes.tick_params(axis='x', colors=p0.get_color())
        ax1.tick_params(axis='y', colors=p1.get_color())
        ax2.tick_params(axis='y', colors=p2.get_color())

        lines = [p0, p1, p2]
        l = axes.legend(lines, [p.get_label() for p in lines])
        for t in l.get_texts():
            t.set_fontsize('xx-small')

        axes.set_xlabel('Time')
        self.formatXticks(axes)
        self.addGrid(axes)

    def ooplot(self):
        axes = self.add_subplot(111)
        color = axes._get_lines.color_cycle
        ax = []
        plots = []
        for i in range(len(self.subfigures) - 1):
            ax.append(axes.twinx())

        p = self.subplot(axes, self.subfigures[0], next(color))
        plots.append(p)
        position = 1.0
        for axes, subfig in zip(ax, self.subfigures[1:]):
            p = self.subplot_twinx(axes, subfig, next(color), position)
            plots.append(p)
            position += 0.1

        l = axes.legend(plots, [p.get_label() for p in plots])
        for t in l.get_texts():
            t.set_fontsize('xx-small')
        axes.set_xlabel('Time')
        self.formatXticks(axes)
        self.addGrid(axes)


    def addLegend(self, axes):
        l = axes.legend(loc=2, handletextpad=0.1, borderpad=0.07,
                        labelspacing=0.07)
        for t in l.get_texts():
            t.set_fontsize('xx-small')

    def formatXticks(self, axes):
        axes.xaxis.set_major_locator(self.dayLocator)
        axes.xaxis.set_minor_locator(self.hourLocator)
        axes.xaxis.set_major_formatter(self.dateFormat)
        axes.format_xdata = self.dateFormat

    def addGrid(self, axes):
        axes.grid(True, which='major', color='k', linestyle='-', linewidth=0.2)

    def addTitle(self, title):
        self.suptitle(title)

def saveFigure(figure, filename):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    canvas = FigureCanvas(figure)
    canvas.print_figure(filename)
