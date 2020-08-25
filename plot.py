import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
#import datetime
#from astropy.time import Time
#from agilepy.utils.AstroUtils import AstroUtils

fermi_binsize = 0.0208333335001091

def time_tt_to_mjd(timett):
        """
        Convert tt to mjd. Tolerance = 0.0000001 s

        Args:
            timett (float): time in tt format

        Returns:
            the input converted in mjd format.
        """
        return (timett / 86400.0) + 53005.0

def time_mjd_to_tt(timemjd):
        """
        Convert mjd to tt. Tolerance = 0.001 s

        Args:
            timemjd (float): time in mjd format

        Returns:
            the input converted in tt format.
        """
        return (timemjd - 53005.0) * 86400.0


def plot(agile_data, fermi_data):
    f = plt.figure()
    ax = f.add_subplot(111)

    #---AGILE----
    t_mjd = []
    tm = (time_tt_to_mjd(agile_data["tstart"]) + time_tt_to_mjd(agile_data["tstop"])) / 2
    yerr = agile_data["rate_error"] * agile_data["exp"]
    tw = tm -  time_tt_to_mjd(agile_data["tstart"])

    ax.errorbar(tm, agile_data["cts"], color="b", label="AGILE", fmt='.', yerr=yerr, xerr=tw)
    
    #---Fermi----
    tstart = fermi_data["Time_MJD"] - fermi_binsize/2
    tstop = tstart + fermi_binsize
    print([tstart, fermi_data["Time_MJD"], tstop])
    fermi_yerr = fermi_data["count_rate_based_error_(cts/s)"] * fermi_data["exposure_(cm^2/s)"]
    ax.errorbar(fermi_data["Time_MJD"], fermi_data["counts"], color="r", label="FERMI", fmt="none", xerr=[fermi_data["Time_MJD"] - tstart,tstop - fermi_data["Time_MJD"]], yerr=fermi_yerr)
    
    
    ax.ticklabel_format(axis="x", useOffset=False)
    ax.set_ylabel('Photon counts')
    ax.set_xlabel("MJD")
    ax.legend(loc='upper right', shadow=True, fontsize='xx-small')



if __name__ == "__main__":

    #--- Parsing args-------
    parser = argparse.ArgumentParser()
    parser.add_argument("--agile", type=str, help="AGILE AP3 Filepath", required=True)
    parser.add_argument("--fermi", type=str, help="FERMI AP3 Filepath", required=True)
    parser.add_argument("--tstart", type=float, help="Tstart in MJD", required=True)
    parser.add_argument("--tstop", type=float, help="Tstop in MJD", required=True)

    args = parser.parse_args()


    #---- Loading data -----
    agile_data = pd.read_csv(args.agile, header=0, sep=" ")
    fermi_data = pd.read_csv(args.fermi, header=0, sep=" ")

    tstart_tt = time_mjd_to_tt(args.tstart)
    tstop_tt = time_mjd_to_tt(args.tstop)
    
    #---- Selecting data
    agile_data = agile_data[agile_data.tstart >= tstart_tt]
    agile_data = agile_data[agile_data.tstop <= tstop_tt]

    fermi_data = fermi_data[fermi_data.Time_MJD >= args.tstart]

    plot(agile_data, fermi_data)

    plt.show()