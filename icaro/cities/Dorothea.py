##############################
##### base_cities.py
##############################

from __future__ import print_function, division

import abc

import tables as tb

from   invisible_cities.cities.base_cities import City

import invisible_cities.core.  fit_functions as fitf
import invisible_cities.reco.  tbl_functions as tbf
import invisible_cities.reco.pmaps_functions as pmp

class FilterCity(City, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def select_event(self, *args):
        pass

    def select_peaks(self, peaks,
                     Emin, Emax,
                     Lmin, Lmax,
                     Hmin, Hmax):
        filtered = copy.deepcopy(peaks)
        for peak_no, (t, E) in peaks.items():
            sel = fitf.in_range(E, Hmin, Hmax)
            if (not Lmin <= np.count_nonzero(sel) < Lmax or
                not Hmin <= np.max(E[sel])        < Hmax or
                not Emin <= np.sum(E[sel])        < Emax):
                del filtered[peak_no]
            else:
                filtered[peak_no] = t[sel], E[sel]
        return filtered

    def select_Si(self, peaks,
                  Nmin, Nmax):
        filtered = copy.deepcopy(peaks)
        for peak_no, sipms in peaks.items():
            if not Nmin <= len(sipms) < Nmax:
                del filtered[peak_no]
        return filtered
        

class DSTCity(City, metaclass=abc.ABCMeta):
    def __init__(self,
                 run_number        = 0       ,
                 files_in          = None    ,
                 file_out          = None    ,
                 compression       = 'ZLIB4' ,
                 nprint            = 10000   ,

                 nwrite            = 10000   ,
                 table_name        = "Events",
                 table_format      = None    ,
                 table_description = ""      ):
        City.__init__(self,
                      run_number  = run_number ,
                      files_in    = files_in   ,
                      file_out    = file_out   ,
                      compression = compression,
                      nprint      = nprint     )

        self.nwrite = nwrite
        self.write  = tb.open_file(self.output_file,
                                   mode = "w",
                                   filters = tbf.filters(self.compression))

        self.table  = self.write.create_table(self.write.root,
                                              table_name,
                                              table_format,
                                              table_description,
                                              tbf.filters(compression))

        self.row    = self.table.row

    @abc.abstractmethod
    def write_evt(self, evt):
        pass

##############################
##### params.py
##############################

# Class instead of namedtuple in order to support default values
class KrEvent:
    def __init__(self):
        self.evt   = -1
        self.T     = -1

        self.nS1   = -1
        self.S1w   = []
        self.S1h   = []
        self.S1e   = []
        self.S1t   = []

        self.nS2   = -1
        self.S2w   = []
        self.S2h   = []
        self.S2e   = []
        self.S2q   = []
        self.S2t   = []

        self.Nsipm = []
        self.DT    = []
        self.Z     = []
        self.X     = []
        self.Y     = []
        self.R     = []
        self.Phi   = []
        self.Xrms  = []
        self.Yrms  = []

##############################
##### nh5.py
##############################

class KrTable(tb.IsDescription):
    event = tb.  Int32Col(pos= 0)
    time  = tb.Float64Col(pos= 1)
    peak  = tb. UInt16Col(pos= 2)
    nS2   = tb. UInt16Col(pos= 3)

    S1w   = tb.Float64Col(pos= 4)
    S1h   = tb.Float64Col(pos= 5)
    S1e   = tb.Float64Col(pos= 6)
    S1t   = tb.Float64Col(pos= 7)

    S2w   = tb.Float64Col(pos= 8)
    S2h   = tb.Float64Col(pos= 9)
    S2e   = tb.Float64Col(pos=10)
    S2q   = tb.Float64Col(pos=11)
    S2t   = tb.Float64Col(pos=12)

    Nsipm = tb.Float64Col(pos=13)
    DT    = tb.Float64Col(pos=14)
    Z     = tb.Float64Col(pos=15)
    X     = tb.Float64Col(pos=16)
    Y     = tb.Float64Col(pos=17)
    R     = tb.Float64Col(pos=18)
    Phi   = tb.Float64Col(pos=19)
    Xrms  = tb.Float64Col(pos=20)
    Yrms  = tb.Float64Col(pos=21)

##############################
##### Dorothea.py
##############################

import sys
import time
import glob
import copy
from   textwrap import dedent

import numpy as np

from invisible_cities.core.system_of_units_c import units
#from invisible_cities.cities.base_cities     import FilterCity
from invisible_cities.core.configure         import configure

#from invisible_cities.reco.nh5    import KrTable
#from invisible_cities.reco.params import KrEvent

class S12Selector(FilterCity):
    def __init__(self,
                 drift_v     = 1 * units.mm/units.mus,

                 S1_Nmin     = 0,
                 S1_Nmax     = 1000,
                 S1_Emin     = 0,
                 S1_Emax     = np.inf,
                 S1_Lmin     = 0,
                 S1_Lmax     = np.inf,
                 S1_Hmin     = 0,
                 S1_Hmax     = np.inf,

                 S2_Nmin     = 0,
                 S2_Nmax     = 1000,
                 S2_Emin     = 0,
                 S2_Emax     = np.inf,
                 S2_Lmin     = 0,
                 S2_Lmax     = np.inf,
                 S2_Hmin     = 0,
                 S2_Hmax     = np.inf,
                 S2_NSIPMmin = 1,
                 S2_NSIPMmax = np.inf):

        FilterCity.__init__(self)

        self.drift_v     = drift_v

        self.S1_Nmin     = S1_Nmin
        self.S1_Nmax     = S1_Nmax
        self.S1_Emin     = S1_Emin
        self.S1_Emax     = S1_Emax
        self.S1_Lmin     = S1_Lmin
        self.S1_Lmax     = S1_Lmax
        self.S1_Hmin     = S1_Hmin
        self.S1_Hmax     = S1_Hmax

        self.S2_Nmin     = S2_Nmin
        self.S2_Nmax     = S2_Nmax
        self.S2_Emin     = S2_Emin
        self.S2_Emax     = S2_Emax
        self.S2_Lmin     = S2_Lmin
        self.S2_Lmax     = S2_Lmax
        self.S2_Hmin     = S2_Hmin
        self.S2_Hmax     = S2_Hmax
        self.S2_NSIPMmin = S2_NSIPMmin
        self.S2_NSIPMmax = S2_NSIPMmax

    def select_S1(self, s1s):
        return self.select_peaks(s1s,
                                 self.S1_Emin, self.S1_Emax,
                                 self.S1_Lmin, self.S1_Lmax,
                                 self.S1_Hmin, self.S1_Hmax)

    def select_S2(self, s2s, sis):
        s2s = self.select_peaks(s2s,
                                self.S2_Emin, self.S2_Emax,
                                self.S2_Lmin, self.S2_Lmax,
                                self.S2_Hmin, self.S2_Hmax)
        sis = self.select_Si(sis,
                             self.S2_NSIPMmin, self.S2_NSIPMmax)

        valid_peaks = set(s2s) & set(sis)
        for d in [s2s, sis]:
            for peak_no in list(d):
                if peak_no not in valid_peaks:
                    del d[peak_no]

        return s2s, sis

    def integrate_charge(self, sipms):
        intq = ( (key, np.sum(value))
                 for key, value in sipms.items())
        return map(np.array, list(zip(*intq)))

    def width(self, t, to_mus=False):
        w = np.max(t) - np.min(t)
        return w * units.ns / units.mus if to_mus else w

    def select_event(self, evt_number, evt_time, S1, S2, Si):
        evt       = KrEvent()
        evt.event = evt_number
        evt.time  = evt_time * 1e-3 # s

        S1     = self.select_S1(S1)
        S2, Si = self.select_S2(S2, Si)

        if (not self.S1_Nmin <= len(S1) < self.S1_Nmax or
            not self.S2_Nmin <= len(S2) < self.S2_Nmax):
            return False, None

        evt.nS1 = len(S1)
        for peak_no, (t, e) in sorted(S1.items()):
            evt.S1w.append(self.width(t))
            evt.S1h.append(np.max(e))
            evt.S1e.append(np.sum(e))
            evt.S1t.append(t[np.argmax(e)])

        evt.nS2 = len(S2)
        for peak_no, (t, e) in sorted(S2.items()):
            evt.S2w.append(self.width(t, to_mus=True))
            evt.S2h.append(np.max(e))
            evt.S2e.append(np.sum(e))
            evt.S2t.append(t[np.argmax(e)])

            IDs, Qs = self.integrate_charge(Si[peak_no])
            xsipms  = self.xs[IDs]
            ysipms  = self.ys[IDs]
            x       = np.average(xsipms, weights=Qs)
            y       = np.average(ysipms, weights=Qs)
            q       = np.sum    (Qs)

            evt.Nsipm.append(len(IDs))
            evt.S2q  .append(q)

            evt.X    .append(x)
            evt.Y    .append(y)

            evt.Xrms .append((np.sum(Qs * (xsipms-x)**2) / (q - 1))**0.5)
            evt.Yrms .append((np.sum(Qs * (ysipms-y)**2) / (q - 1))**0.5)

            evt.R    .append((x**2 + y**2)**0.5)
            evt.Phi  .append(np.arctan2(y, x))

            dt = evt.S2t[-1] - evt.S1t[-1] if len(evt.S2t) > 0 and len(evt.S1t) > 0 else -1

            evt.DT   .append(dt * units.ns / units.mus)
            evt.Z    .append(dt * self.drift_v)

        return True, evt


class KrSelector(S12Selector, DSTCity):
    def __init__(self,
                 run_number  = 0,
                 files_in    = None,
                 file_out    = None,
                 compression = 'ZLIB4',
                 nprint      = 10000,

                 nwrite            = 10000,
                 table_name        = "Events",
                 table_format      = KrTable,
                 table_description = "Selected Kr events",

                 drift_v     = 1 * units.mm/units.mus,

                 S1_Emin     = 0,
                 S1_Emax     = np.inf,
                 S1_Lmin     = 0,
                 S1_Lmax     = np.inf,
                 S1_Hmin     = 0,
                 S1_Hmax     = np.inf,

                 S2_Nmax     = 1,
                 S2_Emin     = 0,
                 S2_Emax     = np.inf,
                 S2_Lmin     = 0,
                 S2_Lmax     = np.inf,
                 S2_Hmin     = 0,
                 S2_Hmax     = np.inf,
                 S2_NSIPMmin = 1,
                 S2_NSIPMmax = np.inf):

        S12Selector.__init__(self,
                             drift_v     = drift_v,

                             S1_Nmin     = 1,
                             S1_Nmax     = 2,
                             S1_Emin     = S1_Emin,
                             S1_Emax     = S1_Emax,
                             S1_Lmin     = S1_Lmin,
                             S1_Lmax     = S1_Lmax,
                             S1_Hmin     = S1_Hmin,
                             S1_Hmax     = S1_Hmax,

                             S2_Nmin     = 1,
                             S2_Nmax     = S2_Nmax,
                             S2_Emin     = S2_Emin,
                             S2_Emax     = S2_Emax,
                             S2_Lmin     = S2_Lmin,
                             S2_Lmax     = S2_Lmax,
                             S2_Hmin     = S2_Hmin,
                             S2_Hmax     = S2_Hmax,
                             S2_NSIPMmin = S2_NSIPMmin,
                             S2_NSIPMmax = S2_NSIPMmax)

        DSTCity.   __init__(self,
                            run_number       ,
                            files_in         ,
                            file_out         ,
                            compression      ,
                            nprint           ,
                            nwrite           ,
                            table_name       ,
                            table_format     ,
                            table_description)

    def write_evt(self, evt):
        for i in range(evt.nS2):
            self.row["event"] = evt.event
            self.row["time" ] = evt.time
            self.row["peak" ] = i
            self.row["nS2"  ] = evt.nS2

            self.row["S1w"  ] = evt.S1w  [0]
            self.row["S1h"  ] = evt.S1h  [0]
            self.row["S1e"  ] = evt.S1e  [0]
            self.row["S1t"  ] = evt.S1t  [0]

            self.row["S2w"  ] = evt.S2w  [i]
            self.row["S2h"  ] = evt.S2h  [i]
            self.row["S2e"  ] = evt.S2e  [i]
            self.row["S2q"  ] = evt.S2q  [i]
            self.row["S2t"  ] = evt.S2t  [i]

            self.row["Nsipm"] = evt.Nsipm[i]
            self.row["DT"   ] = evt.DT   [i]
            self.row["Z"    ] = evt.Z    [i]
            self.row["X"    ] = evt.X    [i]
            self.row["Y"    ] = evt.Y    [i]
            self.row["R"    ] = evt.R    [i]
            self.row["Phi"  ] = evt.Phi  [i]
            self.row["Xrms" ] = evt.Xrms [i]
            self.row["Yrms" ] = evt.Yrms [i]
            self.row.append()

    def run(self, max_evt=-1):
        nevt_in = nevt_out = 0

        exit_file_loop = False
        for filename in self.input_files:
            if nevt_in > self.nprint:
                self.nprint += self.nprint
                print("Opening {}... {} events analyzed".format(filename, nevt_in))
            try:
                S1s, S2s, S2Sis = pmp.load_pmaps(filename)
            except:
                print("##########", filename)
            with tb.open_file(filename, "r") as h5in:
                event_numbers = h5in.root.Run.events.cols.evt_number[:]
                timestamps    = h5in.root.Run.events.cols.timestamp [:]

            # loop
            for evt_number, evt_time in zip(event_numbers, timestamps):
                nevt_in +=1

                S1 = S1s  .get(evt_number, {})
                S2 = S2s  .get(evt_number, {})
                Si = S2Sis.get(evt_number, {})

                ok, evt = self.select_event(evt_number, evt_time,
                                            S1, S2, Si)
                if ok:
                    nevt_out += 1
                    self.write_evt(evt)

                if nevt_in >= max_evt:
                    exit_file_loop = True
                    break

            if exit_file_loop:
                break
        self.table.flush()

        print(dedent("""
                     Number of events in : {}
                     Number of events out: {}
                     Ratio               : {}
                     """.format(nevt_in, nevt_out, nevt_out/nevt_in)))
        return nevt_in, nevt_out, nevt_in/nevt_out

def KRYPTONITE(argv = sys.argv):
    """KsSelector DRIVER"""

    # get parameters dictionary
    CFP = configure(argv)

    #class instance
    kryptonite = KrSelector(run_number  = CFP.RUN_NUMBER,
                            files_in    = sorted(glob.glob(CFP.FILE_IN)),
                            file_out    = CFP.FILE_OUT,
                            compression = CFP.COMPRESSION,
                            nprint      = CFP.NPRINT,
                            nwrite      = CFP.NWRITE,

                            drift_v     = CFP.DRIFT_V,

                            S1_Emin     = CFP.S1_EMIN,
                            S1_Emax     = CFP.S1_EMAX,
                            S1_Lmin     = CFP.S1_LMIN,
                            S1_Lmax     = CFP.S1_LMAX,
                            S1_Hmin     = CFP.S1_HMIN,
                            S1_Hmax     = CFP.S1_HMAX,

                            S2_Nmax     = CFP.S2_NMAX,
                            S2_Emin     = CFP.S2_EMIN,
                            S2_Emax     = CFP.S2_EMAX,
                            S2_Lmin     = CFP.S2_LMIN,
                            S2_Lmax     = CFP.S2_LMAX,
                            S2_Hmin     = CFP.S2_HMIN,
                            S2_Hmax     = CFP.S2_HMAX,
                            S2_NSIPMmin = CFP.S2_NSIPMMIN,
                            S2_NSIPMmax = CFP.S2_NSIPMMAX)

    t0 = time.time()
    nevts = CFP.NEVENTS if not CFP.RUN_ALL else -1
    # run
    nevt_in, nevt_out, ratio = kryptonite.run(max_evt = nevts)
    t1 = time.time()
    dt = t1 - t0

    print("run {} evts in {} s, time/event = {}".format(nevt_in, dt, dt/nevt_in))

    return nevts, nevt_in, nevt_out

if __name__ == "__main__":
    KRYPTONITE(sys.argv)