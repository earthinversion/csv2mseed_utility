#!/home/utpal/miniconda3/envs/dispersion/bin/python
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
import tqdm


info_string = '''
Python utility program to write mseed file from multiple csv files (by Utpal Kumar, IESAS, 2021/05)
'''

PARSER = argparse.ArgumentParser(description=info_string, epilog="csv file format: 'Datetime', 'X', 'Y', 'Z' (2021-04-17 00:00:00.005829,0.00824,-0.01095,1.00362)")


def pbar_desc(pbar, desc):
    pbar.set_description("{:<40}".format(desc))
    

# Conversion factors
GalFactorGlobal = 980.665
days_to_seconds = 60 * 60 * 24

channels = ['BNX', 'BNY', 'BNZ']

chunksize = 10 ** 6
def _read_csv_as_np_array(csvfile, x_arrayInp=None, y_arrayInp=None, z_arrayInp=None, datetime_arrayInp=None, interval=20):    
    x_array = np.array([],dtype=np.float64)
    y_array = np.array([],dtype=np.float64)
    z_array = np.array([],dtype=np.float64)
    datetime_array = np.array([],dtype='datetime64')
    dff = None
    # startDataLength = 0
    if datetime_arrayInp is not None:
        ## takes some data from the previous file
        x_array = np.array(x_arrayInp,dtype=np.float64)
        y_array = np.array(y_arrayInp,dtype=np.float64)
        z_array = np.array(z_arrayInp,dtype=np.float64)
        datetime_array = np.array(datetime_arrayInp,dtype='str')
        dff = pd.DataFrame(columns=['Datetime', "X", "Y", "Z"])
        dff['Datetime'] = datetime_array
        dff['X'] = x_array
        dff['Y'] = y_array
        dff['Z'] = z_array
        
        dff["Datetime"] = pd.to_datetime(dff["Datetime"])
        dff.set_index("Datetime", inplace=True)
        datetime_array = dff.index.values
        # startDataLength = len(datetime_array)
        # print("startDataLength is ",startDataLength)
    count=0
    for df in pd.read_csv(csvfile, chunksize=chunksize, names=[
                        'Datetime', "X", "Y", "Z"], dtype={"Datetime": "str", "X": np.float64,"Y": np.float64,"Z": np.float64}):
        
        df["Datetime"] = pd.to_datetime(df["Datetime"])
        if count==0:
            # print("first line from csv file: ", df.head(1).to_string(header=False,index_names=False))
            starttime = df.loc[0,"Datetime"]
        
        df.set_index("Datetime", inplace=True)
        if dff is not None:
            dataGap = np.abs(int(dff.index.values[-1]-df.index.values[0]))/10**9
            #only use the old data if the gap from the previous data is minimal
            if dataGap <3 * (interval*10**-3):
                df = pd.concat([dff, df])
        df = df.resample(str(interval)+"ms").last().interpolate(method="nearest")

        datetime_array = np.append(datetime_array, df.index.values)
        x_array = np.append(x_array, df['X'].values)
        y_array = np.append(y_array, df['Y'].values)
        z_array = np.append(z_array, df['Z'].values)
        dff = None
        
        count+=1

    # print("last line from csv file: ", df.tail(1).to_string(header=False,index_names=False))

    return (datetime_array[1:], x_array[1:], y_array[1:], z_array[1:]), starttime.to_datetime64()
    
## read the data from the mseed file
def _read_mseed_file(mseedfile):
    sttmp = read(mseedfile)
    data = sttmp[0].data
    starttime = str(sttmp[0].stats.starttime)[:-1]
    endtime = str(sttmp[0].stats.endtime+sttmp[0].stats.delta)[:-1]
    datetimearray = np.arange(starttime,endtime, dtype='datetime64[20ms]')
    print(starttime, len(datetimearray), len(data))


def _write_mseed(network, station, starttime, datetime_array, z_array, x_array, y_array, sample_rate = 50, demean=False, suffix=None,GalFactor=GalFactorGlobal):
    z_array, x_array, y_array = z_array* GalFactor, x_array* GalFactor, y_array* GalFactor #convert from g to Gal

    timediff = datetime_array[-1]-datetime_array[0]
    # print(int(timediff)/10**9)
    # print(np.timedelta64(timediff, 's'))
    num_seconds = int(timediff)/10**9
    # starttimediff = int(starttime-datetime_array[0])/10**9 #nanosecond to second
    stats = {}
    stats['network'] = network
    stats['station'] = station
    # sys.stdout.write(starttime)
    # print('starttime ->',starttime, starttimestr, datetime_array[0], datetime_array[1])
    
    stats['sampling_rate'] = sample_rate
    # datetimeLB = datetime_array[datetime_array>starttime][0]
    starttimestr = np.datetime_as_string(datetime_array[0])
    stats['starttime'] = UTCDateTime(starttimestr)
    npts = int(num_seconds * sample_rate) + 1
    # print(f"datetime_array_end: {datetime_array[-1]}, starttimeOrig: {starttime}, starttimeMseed: {starttimestr}")

    ## resample indexes to use
    if npts<datetime_array.shape[0]:
        random_idxs = np.arange(0, npts+1)
        datetime_array, z_array, x_array, y_array = datetime_array[random_idxs], z_array[random_idxs], x_array[random_idxs], y_array[random_idxs] #resample
    elif npts>datetime_array.shape[0]:
        return
    
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
        if demean:
            stream_to_write[0].detrend("spline", order=3, dspline=500)

        stream_to_write[0].trim(starttime= UTCDateTime(np.datetime_as_string(starttime)))
        if not suffix:
            stream_to_write.write(f"{network}-{station}-{channels[ii]}.mseed", format='MSEED', byteorder='>', reclen=4096)
        else:
            stream_to_write.write(f"{network}-{station}-{channels[ii]}_{suffix}.mseed", format='MSEED', byteorder='>', reclen=4096)
    # print(stream_to_write)
    return (f"{network}-{station}-{channels[0]}.mseed", f"{network}-{station}-{channels[1]}.mseed", f"{network}-{station}-{channels[2]}.mseed")

def main(args):
    x_array = np.array([],dtype=np.float64)
    y_array = np.array([],dtype=np.float64)
    z_array = np.array([],dtype=np.float64)
    datetime_array = np.array([],dtype='datetime64')

    if args.inputs:
        csvfiles = args.inputs
        pbar = tqdm.tqdm(csvfiles)
        starttime0 = None
        try:
            for filecount, csvfile in enumerate(pbar):
                if os.path.exists(csvfile):
                    if csvfile.endswith(".csv"):
                        pbar_desc(pbar, f"Reading file {os.path.basename(csvfile)}")
                        interval =  1/args.sample_rate*1000


                        ## convert to gal or keep in g

                        lengthToLevel = 10
                        if len(datetime_array)>lengthToLevel:
                            datetime_array, x_array, y_array, z_array = datetime_array[-lengthToLevel:], x_array[-lengthToLevel:], y_array[-lengthToLevel:], z_array[-lengthToLevel:]
                            starttime0 = datetime_array[-1]
                            datetime_array = [str(dt) for dt in datetime_array]
                            (datetime_array, x_array, y_array, z_array), starttime =_read_csv_as_np_array(csvfile, interval=interval, datetime_arrayInp=datetime_array, x_arrayInp=x_array, y_arrayInp=y_array, z_arrayInp=z_array)
                        else:
                            (datetime_array, x_array, y_array, z_array), starttime =_read_csv_as_np_array(csvfile, interval=interval)
                        if args.gal:
                            GalFactor = GalFactorGlobal
                        else:
                            GalFactor = 1
                                
                        # pbar.set_description(f"Writing mseed files for {os.path.basename(csvfile)}")
                        # print('first 2',z_array[:2])
                        # print('last 2',z_array[-2:])
                        if starttime0 is not None:
                            starttime1 = starttime0
                        else:
                            starttime1 = starttime

                        # print(starttime, starttime1)
                        mseedX, mseedY, mseedZ = _write_mseed(args.network, args.station, starttime1, datetime_array, z_array, x_array, y_array, sample_rate = args.sample_rate, demean=args.demean, suffix=str(filecount), GalFactor=GalFactor)
                    
                else:
                    sys.stdout.write(f"No file found: {csvfile}\n")
        except Exception as err:
            sys.stdout.write(str(err))





if __name__ == '__main__':
    PARSER.add_argument("-inp",'--inputs', nargs='+', type=str, help="input CSV files to convert to mseed, e.g. network_station_data.csv", required=True)
    PARSER.add_argument("-stn",'--station', type=str, default="XYZ", help="station name, e.g. XYZ")
    PARSER.add_argument("-net", '--network', type=str, default="TW", help="network name, e.g. TW")
    PARSER.add_argument("-sr", '--sample_rate', type=int, default=50, help="sampling rate as integer")
    PARSER.add_argument("-gal", '--gal', type=int, default=1, help="1 for Gal and 0 for g")
    PARSER.add_argument("-dm", '--demean', help="remove mean from the data", action='store_true')

    args = PARSER.parse_args()
    main(args)