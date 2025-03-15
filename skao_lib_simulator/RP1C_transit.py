# docs:
#  - https://casacore.github.io/python-casacore/index.html
#  - https://stellar-h2020.eu/index.php/2023/05/23/introduction-to-lofar-data-processing-tutorial/

# python3 -m venv py3.8
# source py3.8/bin/activate
# pip install numpy matplotlib python-casacore ipython astropy
# ipython
# run RP1C_transit.py

from casacore.tables import table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, AutoDateLocator, AutoDateFormatter, date2num
from astropy.time import Time
from astropy.visualization import time_support


def main(ms_file):
    print(f"Processing {ms_file}")
    MyMS = table(ms_file)

    # only looking at MR01 and MR08, but you can access all baselines
    t18 = MyMS.query('ANTENNA1=1 AND ANTENNA2=8')
    t11 = MyMS.query('ANTENNA1=1 AND ANTENNA2=1')
    t88 = MyMS.query('ANTENNA1=8 AND ANTENNA2=8')

    t_MJDs = t18.getcol('TIME')
    t_MJD = Time(t_MJDs / (3600*24), format='mjd')

    SW = table(ms_file + '::SPECTRAL_WINDOW')
    CHAN_FREQ = SW.getcol('CHAN_FREQ')[0]
    CHAN_WIDTH = SW.getcol('CHAN_WIDTH')[0][0]
    BW = np.array((CHAN_FREQ[0]-CHAN_WIDTH/2, CHAN_FREQ[-1]+CHAN_WIDTH/2))

    for tij in (t18, t11, t88):
        data = tij.getcol('DATA')
        
        fig, (ax_abs, ax_phase) = plt.subplots(ncols=2, nrows=1)
        locator = AutoDateLocator()
        fig.suptitle(f'BW = {BW*1e-6} MHz')
        
        ax_abs.plot(t_MJD.to_datetime(), np.absolute(data[:, :, 0].mean(axis=1)), color='r')
        ax_abs.plot(t_MJD.to_datetime(), np.absolute(data[:, :, 3].mean(axis=1)), color='b')
        ax_abs.set_xlabel('Time')
        ax_abs.set_ylabel('Amp')
        ax_abs.xaxis.set_major_locator(locator)
        ax_abs.xaxis.set_major_formatter( AutoDateFormatter(locator) )
        
        ax_phase.plot(t_MJD.to_datetime(), np.angle(data[:, :, 0]), color='r')
        ax_phase.plot(t_MJD.to_datetime(), np.angle(data[:, :, 3]), color='b')
        ax_phase.set_xlabel('Time')
        ax_phase.set_ylabel('Phase')
        ax_phase.xaxis.set_major_locator(locator)
        ax_phase.xaxis.set_major_formatter( AutoDateFormatter(locator) )
        
        fig.autofmt_xdate()
        
    plt.show()
if __name__ == "__main__":
    main(sys.argv[1])  # Récupère l'argument du fichier MS
