from PyQt5 import QtCore, QtWidgets, uic
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

from PyQt5.uic import loadUiType


qtCreatorFile = "event_maps.ui" # Enter file here.
Ui_MainWindow, QMainWindow = loadUiType(qtCreatorFile)

class QtGui(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        super(QtGui, self).__init__()
        self.setupUi(self)

        self.fig_dict = {}
        self.mplfigs.itemClicked.connect(self.change_figure)

        fig = Figure()
        self.add_mpl(fig)

    def add_figure(self, name, fig):
        test_key = self.fig_dict.pop(name, None)
        self.fig_dict[name] = fig
        if not test_key: # key not in dict
            self.mplfigs.addItem(name)

    def change_figure(self, item):
        text = item.text()
        self.rm_mpl()
        self.add_mpl(self.fig_dict[text])

    def add_mpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,
                self.mplwindow, coordinates=True)
        self.mplvl.addWidget(self.toolbar)

    def rm_mpl(self,):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()
