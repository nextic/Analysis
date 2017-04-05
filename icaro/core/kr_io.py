import tables as tb
import invisible_cities.reco.tbl_functions as tbl
import invisible_cities.reco.nh5 as table_formats
from invisible_cities.reco.pmap_io import event_writer, run_writer,\
     _make_run_event_tables

from invisible_cities.core.system_of_units_c import units
from icaro.core.pmp_functions import s12_features

class KrEventTable(tb.IsDescription):
    """Store for a Kripton event."""
    run       =  tb.  Int32Col(pos=0)
    event     =  tb.  Int32Col(pos=1)
    timestamp =  tb. UInt64Col(pos=2)
    peak      =  tb.  Int16Col(pos=3)
    s1_energy =  tb.Float32Col(pos=4)
    s2_energy =  tb.Float32Col(pos=5)
    s1_epeak  =  tb.Float32Col(pos=6)
    s2_epeak  =  tb.Float32Col(pos=7)
    s1_tpeak  =  tb.Float32Col(pos=8)
    s2_tpeak  =  tb.Float32Col(pos=9)
    s1_width  =  tb.Float32Col(pos=10)
    s2_width  =  tb.Float32Col(pos=11)
    n_s1_pmt  =  tb.  Int32Col(pos=12)
    n_sipm    =  tb.  Int32Col(pos=13)
    Q_tp      =  tb.Float32Col(pos=14)
    x         =  tb.Float32Col(pos=15)
    y         =  tb.Float32Col(pos=16)
    z         =  tb.Float32Col(pos=17)
    r         =  tb.Float32Col(pos=18)
    phi       =  tb.Float32Col(pos=19)


class kr_writer:

    def __init__(self, filename, compression = 'ZLIB4'):
        self._hdf5_file   = tb.open_file(filename, 'w')
        #self._run_tables  = _make_run_event_tables(self._hdf5_file,
        #                                                  compression)
        self._kr_table    = _make_kr_table(self._hdf5_file,
                                           compression)

    def __call__(self, run_number, event_number, peak, timestamp, kr_event):

        kr_event.store(self._kr_table, run_number,
                       event_number, timestamp, peak)

        #run_writer(self._run_tables[0], run_number)
        #event_writer(self._run_tables[1], event_number, timestamp)

    def close(self):
        self._hdf5_file.close()

    @property
    def file(self):
        return self._hdf5_file

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def _make_kr_table(hdf5_file, compression):

    c = tbl.filters(compression)
    kr_group  = hdf5_file.create_group(hdf5_file.root, 'KrEvent')

    MKT = hdf5_file.create_table
    kr_table = MKT(kr_group, 'KrEvent'  , KrEventTable, "KrEvent Table", c)
    kr_table.cols.event.create_index()

    return kr_table


class KrEvent:
    """
    Defines the global features of Kripton Event.
    peak: peak number
    peak goes from tmin to tmax, maximum in tpeak
    width: tmax - tmin
    etot: integrated energy
    emax: energy in tpeak
    er:   emax/etot
    S12F --> {peak:namedtuple('S12Features','tmin tmax tpeak etot epeak')}
    """
    def __init__(self, drift_velocity):
        self.drift_velocity = drift_velocity

    def store(self, table, run_number, event_number, timestamp, peak):
        row = table.row
        row["run"]          = run_number
        row["event"]        = event_number
        row["timestamp"]    = timestamp
        row["peak"]         = peak
        row["s1_energy"]    =  self.s1f.etot
        row["s2_energy"]    =  self.s2f.etot
        row["s1_epeak"]     =  self.s1f.epeak
        row["s2_epeak"]     =  self.s2f.epeak
        row["s1_tpeak"]     =  self.s1f.tpeak
        row["s2_tpeak"]     =  self.s2f.tpeak
        row["s1_width"]     =  self.s1f.width
        row["s2_width"]     =  self.s2f.width
        row["n_s1_pmt"]     =  self.ns1
        row["n_sipm"]       =  self.n_sipm
        row["Q_tp"]         =  self.p.Q
        row["x"]            =  self.p.x
        row["y"]            =  self.p.y
        row["z"]            =  self.p.z
        #row["z"]            =  self.s2f.tpeak - self.s1f.tpeak
        row["r"]            =  self.p.r
        row["phi"]          =  self.p.phi

        row.append()

    def add_s1_features(self, S1, peak_number=0):
        self.s1f = s12_features(S1, peak_number=peak_number)

    def add_ns1_pmt(self, n_match_s1):
        self.ns1 = n_match_s1

    def add_nsipm(self, s2si, peak_number=0):
        self.n_sipm = len(s2si[peak_number].keys())

    def add_s2_features(self, S2, peak_number=0):
        self.s2f = s12_features(S2, peak_number=peak_number)

    def add_coordinates(self, krc):
        self.p =  krc
