"""
Utpal Kumar, 2021/04
Python utility program to write mseed file from csv
"""
import sys, os
from datetime import datetime
import numpy as np
import pandas as pd
import argparse
from obspy.core import UTCDateTime
from obspy import Stream, Trace
import matplotlib.pyplot as plt
from obspy import read
from matplotlib import dates


info_string = '''
Python utility program to write mseed file from csv based on Pandas, Numpy and Obspy (by Utpal Kumar, IESAS, 2021/04)
'''

PARSER = argparse.ArgumentParser(description=info_string, epilog="csv file format: 'Datetime', 'X', 'Y', 'Z' (2021-04-17 00:00:00.005829,0.00824,-0.01095,1.00362)")

# Conversion factors
GalFactor = 980.665
days_to_seconds = 60 * 60 * 24

channels = ['BNX', 'BNY', 'BNZ']

chunksize = 10 ** 6
def _read_csv_as_np_array(csvfile, interval=20, demean=False):

    x_array = np.array([],dtype=np.float64)
    y_array = np.array([],dtype=np.float64)
    z_array = np.array([],dtype=np.float64)
    datetime_array = np.array([],dtype='datetime64')
    count=0
    for df in pd.read_csv(csvfile, chunksize=chunksize, names=[
                        'Datetime', "X", "Y", "Z"], dtype={"Datetime": "str", "X": np.float64,"Y": np.float64,"Z": np.float64}):
        df["Datetime"] = pd.to_datetime(df["Datetime"])
        if count==0:
            starttime = df.loc[0,"Datetime"]
        df.set_index("Datetime", inplace=True)
        df = df.resample(str(interval)+"ms").last().interpolate(method="nearest")

        datetime_array = np.append(datetime_array, df.index.values)
        x_array = np.append(x_array, df['X'].values)
        y_array = np.append(y_array, df['Y'].values)
        z_array = np.append(z_array, df['Z'].values)
        count+=1

    if demean:
        x_array = x_array - np.mean(x_array)
        y_array = y_array - np.mean(y_array)
        z_array = z_array - np.mean(z_array)

    return (datetime_array, x_array, y_array, z_array), starttime.to_datetime64()
    
def _write_mseed(network, station, starttime, datetime_array, z_array, x_array, y_array, sample_rate = 50):
    
    starttimestr = np.datetime_as_string(starttime)
    stats = {}
    stats['network'] = network
    stats['station'] = station
    # sys.stdout.write(starttime)
    stats['starttime'] = UTCDateTime(starttimestr)


    stats['sampling_rate'] = sample_rate
    npts = days_to_seconds * sample_rate

    ## randomly generate indexes to use
    if not npts==datetime_array.shape[0]:
        random_idxs = np.arange(0, npts+1)
        datetime_array, z_array, x_array, y_array = datetime_array[random_idxs], z_array[random_idxs], x_array[random_idxs], y_array[random_idxs] #resample
    
    stats['location'] = '00'
    
    statsZ = stats.copy()
    statsZ['channel'] = 'BNZ'
    
    statsX = stats.copy()
    statsX['channel'] = 'BNX'
    
    statsY = stats.copy()
    statsY['channel'] = 'BNY'
    statslist = [statsX, statsY, statsZ]

    for ii, trace_array in enumerate([x_array, y_array, z_array]):
        trace_to_write = Trace(data=np.array(trace_array*10**6, dtype=np.int32), header=statslist[ii])
                        
        stream_to_write = Stream(traces=[trace_to_write])
        stream_to_write.interpolate(sampling_rate=sample_rate, starttime=stats['starttime'], npts=npts)
        
        stream_to_write.write(f"{network}-{station}-{channels[ii]}.mseed", format='MSEED', byteorder='>', reclen=4096)

def _plot_mseed(network, station, unit="Gal"):
    fig, ax = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
    network = 'TW'
    count = 0
    colors= ['r', 'g', 'b']
    for ii in range(3):
        sttmp = read(f"{network}-{station}-{channels[ii]}.mseed")
        sttmp[0].data = sttmp[0].data/10**6
        ax[ii].plot(sttmp[0].times("matplotlib"), sttmp[0].data, color= colors[ii], lw=0.5, label=f"{network}-{station}-{channels[ii]}")
        ax[ii].set_ylabel(channels[ii]+f" in {unit}")
        count+=1

    hfmt = dates.DateFormatter('%m/%d %H:%M')
    for axx in ax:
        axx.xaxis.set_major_formatter(hfmt)
        axx.grid(True)
        axx.legend(loc=1)
        
    
    ax[2].set_xlabel("Time")
    fig.autofmt_xdate()
    plt.savefig(f'RFidget-plot-{network}-{station}.png', bbox_inches='tight', dpi=300)
    plt.close('all')

def main(args):
    if args.input:
        
        csvfile = os.path.abspath(args.input)
        if os.path.exists(csvfile):
            sys.stdout.write(f"Reading file {csvfile} in chunks...\n")
            try:
                interval =  1/args.sample_rate*1000
                (datetime_array, x_array, y_array, z_array), starttime =_read_csv_as_np_array(csvfile, interval=interval, demean=args.demean)

                ## convert to gal or keep in g
                if args.gal:
                    GalFactor = 980.665
                    outunit = "Gal"
                else:
                    GalFactor = 1
                    outunit = "g"

                z_array, x_array, y_array = z_array* GalFactor, x_array* GalFactor, y_array* GalFactor #convert from g to Gal
                        
                _write_mseed(args.network, args.station, starttime, datetime_array, z_array, x_array, y_array, sample_rate = args.sample_rate)

                sys.stdout.write(f"Finished writing file {args.network}-{args.station}-BNX.mseed\n")
                sys.stdout.write(f"Finished writing file {args.network}-{args.station}-BNY.mseed\n")
                sys.stdout.write(f"Finished writing file {args.network}-{args.station}-BNZ.mseed")
                sys.stdout.write(f"... with sampling rate: {args.sample_rate} Hz\n")

                if args.plot_data:
                    _plot_mseed(args.network, args.station, unit=outunit)
            except Exception as err:
                sys.stdout.write(str(err))

        else:
            sys.stdout.write(f"No file found: {csvfile}")



if __name__ == '__main__':
    PARSER.add_argument("-inp",'--input', type=str, help="input CSV file to convert to mseed, e.g. network_station_data.csv", required=True)
    PARSER.add_argument("-stn",'--station', type=str, default="XYZ", help="station name, e.g. XYZ")
    PARSER.add_argument("-net", '--network', type=str, default="TW", help="network name, e.g. TW")
    PARSER.add_argument("-sr", '--sample_rate', type=int, default=50, help="sampling rate as integer")
    PARSER.add_argument("-gal", '--gal', type=int, default=1, help="1 for Gal and 0 for g")
    PARSER.add_argument("-p", '--plot_data', help="plot the output mseed data", action='store_true')
    PARSER.add_argument("-dm", '--demean', help="remove mean from the data", action='store_true')

    args = PARSER.parse_args()
    main(args)