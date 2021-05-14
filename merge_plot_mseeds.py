"""
Utpal Kumar, 2021/05
Merge and plot multiple mseed data
"""
import os
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from obspy import read
from matplotlib import dates
import tqdm

info_string = '''
Python utility program to merge mseed files and plot the time series and spectrogram (by Utpal Kumar, IESAS, 2021/05)
'''

PARSER = argparse.ArgumentParser(description=info_string, epilog="csv file format: 'Datetime', 'X', 'Y', 'Z' (2021-04-17 00:00:00.005829,0.00824,-0.01095,1.00362)")


def pbar_desc(pbar, desc, update=True):
    pbar.set_description("{:<40}".format(desc))
    if update:
        pbar.update(1)
    pbar.refresh()

def main(args):
    mseedfiles = args.inputs
    # try:
    mseedfile1 = os.path.abspath(mseedfiles[0])
    basepath = os.path.dirname(mseedfile1)
    filepathList = [os.path.abspath(mseed) for mseed in mseedfiles]
    
    network = args.network
    station = args.station

    
    if args.merge_data:
        if args.plot_spectrogram:
            pbar = tqdm.trange(6, desc="Merging MiniSeed Files", leave=True)
        else:
            pbar = tqdm.trange(5, desc="Merging MiniSeed Files", leave=True)


        pbar_desc(pbar, "Reading Files")
        # pbar.refresh()
        for ii, mseedfile in enumerate(filepathList):
            if mseedfile.endswith('.mseed'):
                desc = f"Reading file {os.path.basename(mseedfile)}"
                pbar_desc(pbar, desc, update=False)
                sttmp = read(mseedfile)
                if ii==0:
                    mystream = sttmp
                else:
                    mystream += sttmp
        desc = f"Merging data"
        pbar_desc(pbar, desc)

        # Merge the data together
        mystream.merge(method=1)

        pbar_desc(pbar, "Writing the merged data")
        mystream.write(f"{network}-{station}-merged.mseed", format='MSEED', byteorder='>', reclen=4096)


        desc = "Plotting data"
        pbar_desc(pbar, desc)
        mystream[0].data = mystream[0].data/10**6

        fig, ax = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
        ax.plot(mystream[0].times("matplotlib"), mystream[0].data, color= f"C0", lw=0.5, label=f"{network}-{station}")
        ax.set_ylabel("Data Amplitude")

        hfmt = dates.DateFormatter('%m/%d %H:%M')
        ax.xaxis.set_major_formatter(hfmt)
        ax.grid(True)
        ax.legend(loc=1)
            

        ax.set_xlabel("Time")
        fig.autofmt_xdate()
        plt.savefig(os.path.join(basepath,f'RFidget-plot-{network}-{station}.png'), bbox_inches='tight', dpi=300)
        plt.close('all')
        

        if args.plot_spectrogram:
            from matplotlib.colors import LogNorm
            desc = f"Plotting spectrogram"
            pbar_desc(pbar, desc)

            myArray = mystream[0].data
            streamTimeData = mystream[0].times("matplotlib")
            samplingFrequency = mystream[0].stats.sampling_rate
            fig, ax1 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
            
            powerSpectrum, freqsdata, timedata, im = ax1.specgram(myArray, cmap='jet', Fs=samplingFrequency)
            extent = [timedata.min(), timedata.max(), freqsdata.min(), freqsdata.max()]
            
            mesh = ax1.imshow(powerSpectrum, extent=extent, origin='lower', aspect='auto',cmap='jet', norm=LogNorm(),interpolation='nearest')
            ax1.axis('tight')
            ax1.set(title=f"Power spectrogram {network}-{station}")
            ax1.set_ylabel("Frequency (Hz)")
            ax1.set_xlabel("Time Samples")


            fig.colorbar(mesh, ax=ax1)
            
            plt.savefig(os.path.join(basepath,f'RFidget-spectrogram-{network}-{station}.png'), bbox_inches='tight', dpi=300)
            plt.close('all')
        desc = "Done"
        pbar_desc(pbar, desc)
    else:
        if args.plot_spectrogram:
            pbar = tqdm.trange(len(filepathList)+3, desc="Merging MiniSeed Files", leave=True)
        else:
            pbar = tqdm.trange(len(filepathList)+2, desc="Merging MiniSeed Files", leave=True)
        if args.plot_spectrogram:
            myArray = np.array([],dtype=np.float64)

        # No merge, plot separately
        fig, ax = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
        for ii, mseedfile in enumerate(filepathList):
            pbar_desc(pbar, f"Reading file {os.path.basename(mseedfile)}")
            if mseedfile.endswith('.mseed'):
                sttmp = read(mseedfile)
                # print(sttmp)
                sttmp[0].data = sttmp[0].data/10**6
                if args.plot_spectrogram:
                    myArray = np.append(myArray, sttmp[0].data)
                    samplingFrequency = float(sttmp[0].stats.sampling_rate)
                ax.plot(sttmp[0].times("matplotlib"), sttmp[0].data, color= f"C{ii}", lw=0.5, label=f"{network}-{station}-{ii}")
                ax.set_ylabel("Data Amplitude")
        pbar_desc(pbar, f"Plotting data")

        hfmt = dates.DateFormatter('%m/%d %H:%M')
        ax.xaxis.set_major_formatter(hfmt)
        ax.grid(True)
        ax.legend(loc=1)
            
        ax.set_xlabel("Time")
        fig.autofmt_xdate()
        plt.savefig(os.path.join(basepath,f'RFidget-plot-{network}-{station}.png'), bbox_inches='tight', dpi=300)
        plt.close('all')

        if args.plot_spectrogram:
            from matplotlib.colors import LogNorm
            desc = f"Plotting spectrogram"
            pbar_desc(pbar, desc)
            fig, ax1 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
            
            powerSpectrum, freqsdata, timedata, im = ax1.specgram(myArray, cmap='jet', Fs=samplingFrequency)
            extent = [timedata.min(), timedata.max(), freqsdata.min(), freqsdata.max()]
            
            mesh = ax1.imshow(powerSpectrum, extent=extent, origin='lower', aspect='auto',cmap='jet', norm=LogNorm(),interpolation='nearest')
            ax1.axis('tight')
            ax1.set(title=f"Power spectrogram {network}-{station}")
            ax1.set_ylabel("Frequency (Hz)")
            ax1.set_xlabel("Time Samples")
            fig.colorbar(mesh, ax=ax1)
            plt.savefig(os.path.join(basepath,f'RFidget-spectrogram-{network}-{station}.png'), bbox_inches='tight', dpi=300)
            plt.close('all')
            
        desc = "Done"
        pbar_desc(pbar, desc)




    # except Exception as e:
    #     print(sys.exc_info())


if __name__ == '__main__':
    PARSER.add_argument("-inp",'--inputs', nargs='+', type=str, help="input CSV files to convert to mseed, e.g. network_station_data.csv", required=True)
    PARSER.add_argument("-stn",'--station', type=str, default="XYZ", help="station name, e.g. XYZ")
    PARSER.add_argument("-net", '--network', type=str, default="TW", help="network name, e.g. TW")
    PARSER.add_argument("-m", '--merge_data', help="merge the list of mseed files", action='store_true')
    PARSER.add_argument("-ps", '--plot_spectrogram', help="plot the spectrogram of the input mseed data", action='store_true')


    args = PARSER.parse_args()
    main(args)