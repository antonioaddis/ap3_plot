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


def search_interval(arr1, arr2):
    result = []
    i = 0
    j = 0

    n = len(arr1)
    m = len(arr2)
    while i < n and j < m:
        l = max(arr1[i][0], arr2[j][0])
        r = min(arr1[i][1], arr2[j][1])
        if l < r:
            #print('{', l, ',', r, '}')
            result.append([l,r])
        if arr1[i][1] < arr2[j][1]:
            i += 1
        else:
            j += 1
    return result

def plot(ax, agile_data, fermi_data, arg_lines, plotrate=0):

    #---AGILE----

    tm = (time_tt_to_mjd(agile_data["tstart"]) + time_tt_to_mjd(agile_data["tstop"])) / 2
    #agile_data.loc[agile_data['cts'] == 0, 'rateError'] = 0
    #yerr = agile_data["rateError"]*1e8
    if plotrate == 1:
        yerr = agile_data["rateError"]*1e8
        #print(agile_data["tstart"], agile_data["tstop"], agile_data["cts"], agile_data["exp"], agile_data["rate"]*1e8, agile_data["rateError"]*1e8)
        tw = tm -  time_tt_to_mjd(agile_data["tstart"])
        ax.errorbar(tm, agile_data["rate"]*1e8, color="b", label="AGILE", fmt='.', yerr=yerr, xerr=tw, linewidth=0.8)
        print('AGILE mean', agile_data["rate"].mean()*1e8)
        print('AGILE median', agile_data["rate"].median()*1e8)
        print('AGILE std', agile_data["rate"].std()*1e8)
        agilemean = agile_data["rate"].median()*1e8
        agilestd = agile_data["rate"].std()*1e8
    if plotrate == 0:
        yerr = agile_data["rateError"]*agile_data["exp"]
        #print(agile_data["tstart"], agile_data["tstop"], agile_data["cts"], agile_data["exp"], agile_data["rate"]*1e8, agile_data["rateError"]*1e8)
        tw = tm -  time_tt_to_mjd(agile_data["tstart"])
        ax.errorbar(tm, agile_data["cts"], color="b", label="AGILE", fmt='.', yerr=yerr, xerr=tw, linewidth=0.8)
        print('AGILE mean', agile_data["cts"].mean())
        print('AGILE median', agile_data["cts"].median())
        print('AGILE std', agile_data["cts"].std())
        agilemean = agile_data["cts"].median()
        agilestd = agile_data["cts"].std()

    ax.axhline(agilemean, linestyle='solid', color='b', linewidth=0.5)
    ax.axhline(agilemean + 1 * agilestd, linestyle='dotted', color='b', linewidth=0.5)
    ax.axhline(agilemean + 2 * agilestd, linestyle='dashed', color='b', linewidth=0.5)
    ax.axhline(agilemean + 3 * agilestd, linestyle='dashdot', color='b', linewidth=1)

    #---Fermi----
    # tstart = fermi_data["Time_MJD"] - fermi_binsize/2
    # tstop = tstart + fermi_binsize
    # #print([tstart, fermi_data["Time_MJD"], tstop])
    # fermi_yerr = fermi_data["count_rate_based_error_(cts/s)"] * fermi_data["exposure_(cm^2/s)"]
    tmFermi = (time_tt_to_mjd(fermi_data["tstart"]) + time_tt_to_mjd(fermi_data["tstop"])) / 2
    #fermi_data.loc[fermi_data['cts'] == 0, 'rateError'] = 0
    if plotrate == 1:
        fermi_data.loc[fermi_data['rateError'] > 1000e-08, 'rateError'] = 0
        fermi_data.loc[fermi_data['rateError'] > 1000e-08, 'rate'] = 0
        fermi_data.loc[fermi_data['rate'] > 4000e-08, 'rate'] = 0
        yerrFermi = fermi_data["rateError"]*1e8
        #print(fermi_data["tstart"], fermi_data["tstop"], fermi_data["rateError"], fermi_data["rate"])
        twFermi = tmFermi -  time_tt_to_mjd(fermi_data["tstart"])
        #ax.errorbar(fermi_data["Time_MJD"], fermi_data["count_rate_(cts/s)"]*1e8, color="r", label="FERMI", fmt="none", xerr=[fermi_data["Time_MJD"] - tstart,tstop - fermi_data["Time_MJD"]], yerr=fermi_yerr)
        ax.errorbar(tmFermi, fermi_data["rate"]*1e8, color="r", label="FERMI", fmt="none", yerr=yerrFermi, xerr=twFermi, linewidth=0.8)
        print('Fermi mean', fermi_data["rate"].mean()*1e8)
        print('Fermi median', fermi_data["rate"].median()*1e8)
        print('Fermi std', fermi_data["rate"].std()*1e8)
        fermimean = fermi_data["rate"].median()*1e8
        fermistd = fermi_data["rate"].std()*1e8
    if plotrate == 0:
        #fermi_data.loc[fermi_data['rateError'] > 1000e-08, 'rateError'] = 0
        #fermi_data.loc[fermi_data['rateError'] > 1000e-08, 'rate'] = 0
        yerrFermi = fermi_data["rateError"]*fermi_data["exp"]
        #print(fermi_data["tstart"], fermi_data["tstop"], fermi_data["rateError"], fermi_data["rate"])
        twFermi = tmFermi -  time_tt_to_mjd(fermi_data["tstart"])
        #ax.errorbar(fermi_data["Time_MJD"], fermi_data["count_rate_(cts/s)"]*1e8, color="r", label="FERMI", fmt="none", xerr=[fermi_data["Time_MJD"] - tstart,tstop - fermi_data["Time_MJD"]], yerr=fermi_yerr)
        ax.errorbar(tmFermi, fermi_data["cts"], color="r", label="FERMI", fmt="none", yerr=yerrFermi, xerr=twFermi, linewidth=0.8)
        print('Fermi mean', fermi_data["cts"].mean())
        print('Fermi median', fermi_data["cts"].median())
        print('Fermi std', fermi_data["cts"].std())
        fermimean = fermi_data["cts"].mean()
        fermistd = fermi_data["cts"].std()

    ax.axhline(fermimean, linestyle='solid', color='r', linewidth=0.5)
    ax.axhline(fermimean + 1 * fermistd, linestyle='dotted', color='r', linewidth=0.5)
    ax.axhline(fermimean + 2 * fermistd, linestyle='dashed', color='r', linewidth=0.5)
    ax.axhline(fermimean + 3 * fermistd, linestyle='dashdot', color='r', linewidth=1)

    time_diff = fermi_data["tstop"] - fermi_data["tstart"]
    print("Total time in GTI(bottom plot)" ,time_diff.sum())
        



    ax.ticklabel_format(axis="x", useOffset=False)
    if plotrate == 0:
        ax.set_ylabel('Photon counts')
    else:
       ax.set_ylabel('Rate') 
    ax.set_xlabel("MJD")
    ax.legend(loc='upper right', shadow=True, fontsize='xx-small')

def plot_offaxis(ax1, ax2, path, tstart, tstop, zmax, step, t0, arg_lines):

    try:
        agl_meantime, agl_separation = np.loadtxt(path+'/time_vs_separation_agile.txt', unpack=True)
    except:
        return

    agl_filt = agl_meantime[(agl_meantime > tstart) & (agl_meantime < tstop)]
    agl_sep_filt = agl_separation[(agl_meantime > tstart) & (agl_meantime < tstop)]


    lat_meantime, lat_separation = np.loadtxt(path+'/time_vs_separation_fermi.txt', unpack=True)

    lat_filt = lat_meantime[(lat_meantime > tstart) & (lat_meantime < tstop)]
    lat_sep_filt = lat_separation[(lat_meantime > tstart) & (lat_meantime < tstop)]

    ax1.plot(agl_filt - t0, agl_sep_filt, color='blue', label='AGILE', linewidth=0.5)

    ax1.plot(lat_filt - t0, lat_sep_filt, color='red', label='Fermi', linewidth=1.5)

    #-----green boxes----
    lat_filt2 = []
    for i in range(len(lat_filt) - 1):

        if (lat_filt[i+1] - lat_filt[i]) * 86400 >= 300:

            print("Green box in: ", lat_filt[i], lat_filt[i+1])

            lat_filt2.append([lat_filt[i], lat_filt[i+1]])

            ax1.axvline(lat_filt[i], linestyle='--', color='g', linewidth=0.5)
            ax1.axvline(lat_filt[i+1], linestyle='--', color='g', linewidth=0.5)
    ######


    #######------GTI------###
    found = False
    total_s_in_gti = 0
    gti_list = []
    for l, s in zip(lat_filt, lat_sep_filt):

        if not found and s <= zmax:
            found = True
            gti_time = l * 86400
            #ax.axvline(l, linestyle='--', color='k', linewidth=0.5)
            l1 = l

        if found and s >= zmax:
            found = False
            gti_time = (l*86400) - gti_time
            total_s_in_gti += gti_time
            #ax.axvline(l, linestyle='--', color='k', linewidth=0.5)
            l2 = l
            gti_list.append([l1,l2])
    ######

    
    
    #####-----cleaning------
    result = search_interval(lat_filt2, gti_list)

    for l in result:
        #ax.axvline(l[0], linestyle='--', color='k', linewidth=1)
        #ax.axvline(l[1], linestyle='--', color='k', linewidth=1)

        seconds = (l[1] - l[0]) * 86400

        total_s_in_gti = total_s_in_gti - seconds

    for lines in gti_list:
        ax1.axvspan(xmin=lines[0], xmax=lines[1], facecolor='k', alpha=0.1)
        ax2.axvspan(xmin=lines[0], xmax=lines[1], facecolor='k', alpha=0.1) #bottom plot

    for lines in result:
        ax1.axvspan(xmin=lines[0], xmax=lines[1], facecolor='white')
        ax2.axvspan(xmin=lines[0], xmax=lines[1], facecolor='white') #bottom plot

    ######

    print("Total time in GTI", total_s_in_gti)

    try:
        for i in range(0,len(arg_lines),2):

            #ax1.axvline(arg_lines[i], linestyle='--', color='k', linewidth=0.5)
            #ax1.axvline(arg_lines[i+1], linestyle='--', color='k', linewidth=0.5)
            #ax2.axvline(arg_lines[i], linestyle='--', color='k', linewidth=0.5)
            #ax2.axvline(arg_lines[i+1], linestyle='--', color='k', linewidth=0.5)


            ax1.axvspan(xmin=arg_lines[i], xmax=arg_lines[i+1], facecolor='y', alpha=0.1)
            ax2.axvspan(xmin=arg_lines[i], xmax=arg_lines[i+1], facecolor='y', alpha=0.1)
    except:
        print("No lines")
    

    ax1.set_ylim(0., zmax+5.0)
    #ax.set_xlim((tstart - t0)-0.2, (tstop-t0)+0.2)
    ax1.set_xlabel('MJD')

    ax1.ticklabel_format(axis="x", useOffset=False)

    ax1.legend(loc='lower right', shadow=True, fontsize='xx-small')

    ax1.set_ylabel('off-axis angle [$^{\\circ}$]')

    #ax.set_xlim(np.min(agl_filt-t0), np.max(agl_filt-t0))
    ax1.set_title(str(zmax)+'_'+str(tstart)+'_'+str(tstop))

def checkSignificance(fermi, tstart, tstop):
    fermi_data = pd.read_csv(fermi, header=0, sep=" ")
    ntrials = 0
    nsig = 0

    for time in range(int(tstart), int(tstop)):
        #print(time)
        tstart_tt = time_mjd_to_tt(time)
        tstop_tt = time_mjd_to_tt(time+1)
        fermi_data2 = fermi_data[fermi_data.tstart >= tstart_tt]
        fermi_data2 = fermi_data2[fermi_data.tstop <= tstop_tt]
        fermimean = fermi_data2["cts"].mean()
        fermistd = fermi_data2["cts"].std()

        n=0
        for cts in fermi_data2["cts"]:
            ntrials = ntrials + 1
            #print(time, time+1, cts, fermimean, fermistd, fermimean + 3 * fermistd, cts >= (fermimean + 3 * fermistd))
            if cts >= (fermimean + 5 * fermistd):
                print("####")
                print(fermi_data2["tstart"])
                #print(fermi_data2["tstart"][n])
                nsig = nsig + 1
                break
            n = n + 1

        

    print("ntrials ", ntrials)
    print("nsig ", nsig)

def main(agile, fermi, tstart, tstop, path, lines, plotrate=0):
    #---- Loading data -----
    agile_data = pd.read_csv(agile, header=0, sep=" ")
    fermi_data = pd.read_csv(fermi, header=0, sep=" ")
    #print(agile_data)
    tstart_tt = time_mjd_to_tt(tstart)
    tstop_tt = time_mjd_to_tt(tstop)
    #print(tstart_tt, tstop_tt)

    #---- Selecting data
    agile_data = agile_data[agile_data.tstart >= tstart_tt]
    agile_data = agile_data[agile_data.tstop <= tstop_tt]
    #print(agile_data["rateError"])
    fermi_data = fermi_data[fermi_data.tstart >= tstart_tt]
    fermi_data = fermi_data[fermi_data.tstop <= tstop_tt]

    #------Plotting data
    f, (ax1, ax2) = plt.subplots(2, figsize=(12.18,10))
    plot_offaxis(ax1, ax2, path, tstart, tstop, 60, 1, 0, lines)
    plot(ax2, agile_data, fermi_data, lines, plotrate)

    plt.show()
    f.savefig('merged_plot_'+str(tstart)+'_'+str(tstop)+'.'+str('pdf'), format="pdf")

if __name__ == "__main__":
    #checkSignificance("~/projects/IGRJ17354-3255/2020/LS/FERMI/FER2_TBS300_R2_BARC0_EB100-10000.ap.ap4", 54732, 54733)
    #--- Parsing args-------
    parser = argparse.ArgumentParser()
    parser.add_argument("--agile", type=str, help="AGILE AP3 Filepath", required=True)
    parser.add_argument("--fermi", type=str, help="FERMI AP3 Filepath", required=True)
    parser.add_argument("--tstart", type=float, help="Tstart in MJD", required=True)
    parser.add_argument("--tstop", type=float, help="Tstop in MJD", required=True)
    parser.add_argument("--path_offaxis", type=str, help="offaxis directory path", required=True)
    parser.add_argument("--lines", type=float, nargs="+", help="vertical lines", required=False)

    args = parser.parse_args()

    #main(args.agile, args.fermi, args.tstart, args.tstop, args.path_offaxis, args.lines, plotrate=0)
   
