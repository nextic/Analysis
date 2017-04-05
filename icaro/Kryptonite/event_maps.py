import sys
import numpy as np

from PyQt5 import QtCore, QtWidgets, uic
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

from PyQt5.uic import loadUiType


qtCreatorFile = "event_maps.ui" # Enter file here.
Ui_MainWindow, QMainWindow = loadUiType(qtCreatorFile)

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        super(Main, self).__init__()
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



if __name__ == '__main__':
    import sys
    #from PyQt5 import QtGui
    import numpy as np

    fig1 = Figure()
    ax1f1 = fig1.add_subplot(111)
    ax1f1.plot(np.random.rand(5))

    fig2 = Figure()
    ax1f2 = fig2.add_subplot(121)
    ax1f2.plot(np.random.rand(5))
    ax2f2 = fig2.add_subplot(122)
    ax2f2.plot(np.random.rand(10))

    fig3 = Figure()
    ax1f3 = fig3.add_subplot(111)
    ax1f3.pcolormesh(np.random.rand(20,20))

    app = QtWidgets.QApplication(sys.argv)
    #app = QtGui.QApplication(sys.argv)
    main = Main()
    main.add_figure('One plot', fig1)
    main.add_figure('Two plots', fig2)
    main.add_figure('Pcolormesh', fig3)
    main.show()
    sys.exit(app.exec_())
def event_pmaps(P,
                pmtrwf, sipmrwf,
                s1par, s2par, thr,
                s12f,
                deconv_n_baseline  = 48000,
                deconv_thr_trigger =     5,
                s1_match_pmt       = 1 * units.mus,
                event=0):
    """Compute Event pmaps and pmaps features.
    input:
    P           : CalibVectors named tuple.
                  ('CalibVectors',
                  'channel_id coeff_blr coeff_c adc_to_pes pmt_active')
    pmtrwf      : raw waveform for pmts
    sipmrwf     : raw waveform for SiPMs
    thr         : threshold for PMAP searches.
                  ('ThresholdParams',
                  'thr_s1 thr_s2 thr_MAU thr_sipm thr_SIPM')
    s12f        : named tuple containing features of S1 and S2
                  S12Features = namedtuple('S12Features',
                    's1f s1_pmt_f s1pf s2f ')
    event       : event number
    """

    RWF = pmtrwf[event]
    CWF                  = blr.deconv_pmt(RWF,
                                          P.coeff_c,
                                          P.coeff_blr,
                                          n_baseline  = deconv_n_baseline,
                                          thr_trigger = deconv_thr_trigger)

    # calibrated sum

    csum, csum_mau       = cpf.calibrated_pmt_sum(CWF,
                                                  P.adc_to_pes,
                                                  pmt_active   = P.pmt_active,
                                                  n_MAU        = 100,
                                                  thr_MAU      = thr.thr_MAU)


    CAL_PMT, CAL_PMT_MAU = cpf.calibrated_pmt_mau(CWF,
                                                  P.adc_to_pes,
                                                  pmt_active = P.pmt_active,
                                                  n_MAU      = 100,
                                                  thr_MAU    = thr.thr_MAU)


    s1_ene, s1_indx      = cpf.wfzs(              csum_mau,
                                                  threshold  =thr.thr_s1)
    S1                   =  cpf.find_S12(         s1_ene,
                                                  s1_indx,
                                                  **s1par._asdict())

    print_s12(S1)

    if len(S1) == 0:
        return EventPmaps.not_s1

    s2_ene, s2_indx = cpf.wfzs(csum, threshold=thr.thr_s1)
    S2    = cpf.find_S12(s2_ene, s2_indx, **s2par._asdict())

    print_s12(S2)
    if len(S2) == 0:
        return EventPmaps.not_s2

    for i in S1:
        print('S1: adding features for peak number {}'.format(i))
        s12f.s1f.add_features(self, S1, peak_number=i)

        # Search S1 for PMTs where we have found S1 for sum
        t = S1[i][0]
        tmin = t[0] - s1_match_pmt
        tmax = t[-1] + s1_match_pmt

        s1par_PMT = S12Params(tmin=tmin, tmax=tmax, lmin=3, lmax=20,
                              stride=4, rebin=False)

        PMT_S1 = {}
        for pmt in P.pmt_active:
            s1_ene, s1_indx = cpf.wfzs(CAL_PMT_MAU[pmt], threshold=0.1)
            PMT_S1[pmt] = cpf.find_S12(s1_ene, s1_indx, **s1par_PMT._asdict())
            s12f.s1_pmt_f.add_features(self, PMT_S1[pmt], peak_number=i)


        nm = compare_S1(S1, PMT_S1, peak=i)
        print('number of S1 - S1_PMT matches = {}'.format(nm))
        s12f.s1f.add_number_of_matches(nm)

    for i in S2:
        print('S2: adding features for peak number {}'.format(i))
        s12f.s2.add_features(S2, peak_number=i)

    # S1p
    t = S2[0][0]
    tmin = t[-1] + 5 * units.mus
    print('S1p search starts at {}'.format(tmin))

    s1p_params = S12Params(tmin=tmin, tmax=1300*units.mus, stride=4, lmin=4,
                           lmax=20, rebin=False)
    s1_ene, s1_indx = cpf.wfzs(csum_mau, threshold=0.1)
    S1p =  cpf.find_S12(s1_ene, s1_indx, **s1p_params._asdict())
    for i in S1p:
        s12f.s1pf.add_features(S1p, peak_number=i)


    if len(S1) != 1:
        return EventPmaps.s1_not_1

    if len(S2) != 1:
        return EventPmaps.s2_not_1  # TODO: treat the case with more than 1 s2


    dt = s2f.tpeak - s1f.tpeak

    print('drif time = {} mus'.format(dt/units.mus))

    print('***S2Si****')
    sipm = cpf.signal_sipm(sipmrwf[event], adc_to_pes_sipm, thr=3.5*units.pes, n_MAU=100)
    SIPM = cpf.select_sipm(sipm)
    S2Si = pf.sipm_s2_dict(SIPM, S2, thr=20*units.pes)

    return (CalibratedSum(csum=csum, csum_mau=csum_mau),
            PMaps(S1=S1, S2=S2, S2Si=S2Si))
