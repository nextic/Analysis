"""
Compute PMAPS and PMAP features
"""
import sys
import numpy as np

from   invisible_cities.database import load_db
import invisible_cities.sierpe.blr as blr
import invisible_cities.reco.peak_functions_c as cpf
import invisible_cities.reco.peak_functions as pf
from   invisible_cities.core.system_of_units_c import units
from   invisible_cities.core.core_functions import loc_elem_1d
from   invisible_cities.reco.params import S12Params, ThresholdParams,\
                CalibratedSum, PMaps, CalibVectors, DeconvParams

from icaro.core.pmp_functions import s12_features, print_s12, print_s2si, compare_S1


class EventPmaps:
    """Compute event pmaps

    calib_vectors : named tuple.
                    ('CalibVectors',
                    'channel_id coeff_blr coeff_c adc_to_pes pmt_active')

    deconv_params : named tuple.
                    ('DeconvParams', 'n_baseline thr_trigger')

    """

    def __init__(self, run_number,
                 s1par, s2par, thr,
                 verbose=True):
        """
        input:

        pmtrwf        : raw waveform for pmts
        sipmrwf       : raw waveform for SiPMs
        s1par, s2par  : named tuples
                            ('S12Params' , 'tmin tmax stride lmin lmax rebin')
        thr           : named tuple.
                          ('ThresholdParams',
                          'thr_s1 thr_s2 thr_MAU thr_sipm thr_SIPM')
        verbose       : to make it talk.

        """

        self.run_number = run_number
        DataPMT = load_db.DataPMT(run_number)
        DataSiPM = load_db.DataSiPM(run_number)
        self.xs = DataSiPM.X.values
        self.ys = DataSiPM.Y.values
        pmt_active = np.nonzero(DataPMT.Active.values)[0].tolist()

        self.P = CalibVectors(channel_id = DataPMT.ChannelID.values,
                              coeff_blr = abs(DataPMT.coeff_blr   .values),
                              coeff_c = abs(DataPMT.coeff_c   .values),
                              adc_to_pes = abs(DataPMT.adc_to_pes.values),
                              adc_to_pes_sipm = abs(DataSiPM.adc_to_pes.values),
                              pmt_active = pmt_active)

        self.D             = DeconvParams(n_baseline = 48000,
                                          thr_trigger = 5)
        self.s1par         = s1par
        self.s2par         = s2par
        self.thr           = thr

        self.verbose = verbose

    def calibrated_pmt_and_csum(self, event, pmtrwf):
        """Calibrated PMTs and sum."""
        self.RWF = pmtrwf[event]

        self.CWF = blr.deconv_pmt(self.RWF,
                             self.P.coeff_c,
                             self.P.coeff_blr,
                             self.P.pmt_active,
                             n_baseline  = self.D.n_baseline,
                             thr_trigger = self.D.thr_trigger)


        self.csum, self.csum_mau = cpf.calibrated_pmt_sum(self.CWF,
                                            self.P.adc_to_pes,
                                            pmt_active = self.P.pmt_active,
                                            n_MAU      = 100,
                                            thr_MAU    = self.thr.thr_MAU)

        self.CAL_PMT, self.CAL_PMT_MAU  =  cpf.calibrated_pmt_mau(self.CWF,
                                            self.P.adc_to_pes,
                                            pmt_active = self.P.pmt_active,
                                            n_MAU      = 100,
                                            thr_MAU    = self.thr.thr_MAU)

        return np.sum(self.csum)

    def calibrated_sipm(self, event, sipmrwf):
        """Calibrated SiPMs"""

        sipm = cpf.signal_sipm(sipmrwf[event], self.P.adc_to_pes_sipm,
                               thr=self.thr.thr_sipm, n_MAU=100)
        self.SIPM = cpf.select_sipm(sipm)
        return len(self.SIPM)

    def find_s1(self):
        """Compute S1."""
        s1_ene, s1_indx = cpf.wfzs(self.csum_mau, threshold  =self.thr.thr_s1)
        self.S1         = cpf.find_S12(s1_ene, s1_indx, **self.s1par._asdict())
        self.s1f = {}
        if self.verbose:
            print_s12(self.S1)

        for peak in self.S1:
            self.s1f[peak] = s12_features(self.S1, peak_number=peak)

        return len(self.S1)

    def find_s2(self):
        """Compute S2."""
        s2_ene, s2_indx = cpf.wfzs(self.csum, threshold=self.thr.thr_s2)
        self.S2         = cpf.find_S12(s2_ene, s2_indx, **self.s2par._asdict())
        self.s2f = {}
        if self.verbose:
            print_s12(self.S2)

        for peak in self.S2:
            self.s2f[peak] = s12_features(self.S2, peak_number=peak)

        return len(self.S2)

    def find_ns1(self, s1_delta=100*units.ns):
        """Compute matches between csum S1 and PMT S1 for peak=0"""

        s1par_PMT = S12Params(tmin=self.s1f[0].tmin - s1_delta,
                              tmax=self.s1f[0].tmax + s1_delta,
                              lmin=3, lmax=30, stride=4, rebin=False)

        PMT_S1 = {}

        for pmt in self.P.pmt_active:
            s1_ene, s1_indx = cpf.wfzs(self.CAL_PMT_MAU[pmt], threshold=0.1)
            PMT_S1[pmt] = cpf.find_S12(s1_ene, s1_indx, **s1par_PMT._asdict())


        n_match_s1 = compare_S1(self.S1, PMT_S1, peak=0, tol=0.25*units.mus)

        return n_match_s1

    def find_s2si(self):
        """Compute S2Si"""

        self.S2Si = pf.sipm_s2_dict(self.SIPM, self.S2, thr=self.thr.thr_SIPM)
        if self.verbose:
            print_s2si(self.S2Si)

        return len(self.S2Si)
