import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

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

def plot(ax, agile_data, fermi_data, line1, line2):

    #---AGILE----
    t_mjd = []
    tm = (time_tt_to_mjd(agile_data["tstart"]) + time_tt_to_mjd(agile_data["tstop"])) / 2
    yerr = agile_data["rate_error"] * agile_data["exp"]
    print(agile_data["tstart"], agile_data["tstop"], agile_data["rate_error"], agile_data["rate"])
    tw = tm -  time_tt_to_mjd(agile_data["tstart"])

    ax.errorbar(tm, agile_data["rate"]*1e8, color="b", label="AGILE", fmt='.', yerr=yerr, xerr=tw)
    
    #---Fermi----
    tstart = fermi_data["Time_MJD"] - fermi_binsize/2
    tstop = tstart + fermi_binsize
    #print([tstart, fermi_data["Time_MJD"], tstop])
    fermi_yerr = fermi_data["count_rate_based_error_(cts/s)"] * fermi_data["exposure_(cm^2/s)"]
    
    ax.errorbar(fermi_data["Time_MJD"], fermi_data["count_rate_(cts/s)"]*1e8, color="r", label="FERMI", fmt="none", xerr=[fermi_data["Time_MJD"] - tstart,tstop - fermi_data["Time_MJD"]], yerr=fermi_yerr)
    
    if line1 and line2:
        ax.axvline(line1, linestyle='--', color='k', linewidth=0.5)
        ax.axvline(line2, linestyle='--', color='k', linewidth=0.5)
    
    ax.ticklabel_format(axis="x", useOffset=False)
    ax.set_ylabel('Photon counts')
    ax.set_xlabel("MJD")
    ax.legend(loc='upper right', shadow=True, fontsize='xx-small')

def plot_offaxis(ax, path, tstart, tstop, zmax, step, t0, line1, line2):

    agl_meantime, agl_separation = np.loadtxt(path+'/time_vs_separation_agile.txt', unpack=True)
    
    agl_filt = agl_meantime[(agl_meantime > tstart) & (agl_meantime < tstop)]
    agl_sep_filt = agl_separation[(agl_meantime > tstart) & (agl_meantime < tstop)]


    lat_meantime, lat_separation = np.loadtxt(path+'/time_vs_separation_fermi.txt', unpack=True)

    lat_filt = lat_meantime[(lat_meantime > tstart) & (lat_meantime < tstop)]
    lat_sep_filt = lat_separation[(lat_meantime > tstart) & (lat_meantime < tstop)]

    ax.plot(agl_filt - t0, agl_sep_filt, color='blue', label='AGILE')

    ax.plot(lat_filt - t0, lat_sep_filt, color='red', label='Fermi')
        
    ax.axvline(line2, linestyle='--', color='k', linewidth=0.5)


    ax.set_ylim(0., zmax+5.0)
    #ax.set_xlim((tstart - t0)-0.2, (tstop-t0)+0.2)
    ax.set_xlabel('MJD')

    ax.ticklabel_format(axis="x", useOffset=False)

    ax.legend(loc='lower right', shadow=True, fontsize='xx-small')

    ax.set_ylabel('off-axis angle [$^{\\circ}$]')

    #ax.set_xlim(np.min(agl_filt-t0), np.max(agl_filt-t0))
    ax.set_title(str(zmax)+'_'+str(tstart)+'_'+str(tstop))


if __name__ == "__main__":

    #--- Parsing args-------
    parser = argparse.ArgumentParser()
    parser.add_argument("--agile", type=str, help="AGILE AP3 Filepath", required=True)
    parser.add_argument("--fermi", type=str, help="FERMI AP3 Filepath", required=True)
    parser.add_argument("--tstart", type=float, help="Tstart in MJD", required=True)
    parser.add_argument("--tstop", type=float, help="Tstop in MJD", required=True)
    parser.add_argument("--path_offaxis", type=str, help="offaxis directory path", required=True)
    parser.add_argument("--line1", type=float, help="vertical line1", required=False)
    parser.add_argument("--line2", type=float, help="vertical line2", required=False)

    args = parser.parse_args()


    #---- Loading data -----
    agile_data = pd.read_csv(args.agile, header=0, sep=" ")
    fermi_data = pd.read_csv(args.fermi, header=0, sep=" ")
    print(agile_data)
    tstart_tt = time_mjd_to_tt(args.tstart)
    tstop_tt = time_mjd_to_tt(args.tstop)
    print(tstart_tt, tstop_tt)
    
    #---- Selecting data
    agile_data = agile_data[agile_data.tstart >= tstart_tt]
    agile_data = agile_data[agile_data.tstop <= tstop_tt]
    print(agile_data["rate"]*1e8)
    fermi_data = fermi_data[fermi_data.Time_MJD >= args.tstart]
    fermi_data = fermi_data[fermi_data.Time_MJD <= args.tstop]

    f, (ax1, ax2) = plt.subplots(2)

    plot_offaxis(ax1, args.path_offaxis, args.tstart, args.tstop, 60, 1, 0, args.line1, args.line2)
    plot(ax2, agile_data, fermi_data, args.line1, args.line2)

    plt.show()