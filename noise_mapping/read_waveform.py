#!/usr/bin/env python

#RMS 2018
#Read binary files and metadata
#We wil use this script to to processing on the noise data in each grid cell in our model

import struct 
import pickle
import pandas as pd
import numpy as np
import glob
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from datetime import datetime
from scipy.signal import butter, lfilter, freqz
import time


#Functions for lowpass filtering
#-----------------------------------------
def butter_lowpass(cutoff, fs, order=2):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=2):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
#-----------------------------------------


def resample(X,Y,Z,T,Tnew,resamp=25):

    '''
    Take in X,Y,Z and Y vectors from smartphone and resample
    along the vector Tnew. The results will have the same start time, but
    will be have a regular sample rate given by resamp
    '''

    #lowpass filter to avoid artifacts
    sample_rate = round(1.0/np.median(np.diff(T)))
    X = butter_lowpass_filter(X,20.0,sample_rate)
    Y = butter_lowpass_filter(Y,20.0,sample_rate)
    Z = butter_lowpass_filter(Z,20.0,sample_rate)

    #Generate interpolation functions for the traces
    fx = interp1d(T,X)
    fy = interp1d(T,Y)
    fz = interp1d(T,Z)

    #Interpolate the data to new sample rate
    Xnew = fx(Tnew)
    Ynew = fy(Tnew)
    Znew = fz(Tnew)

    #integrate to velocity
    Xnew = 
 
    return Xnew,Ynew,Znew

def main():

    binary_data_files = glob.glob('data_cell*.bin')
    metadata_files = glob.glob('metadata_cell*.pkl')

    #conversion from g to m/s^2
    g2ms=9.80655

    nfreq = 25
    nsec = 300
    Tnew = np.arange(0,nsec+(1.0/nfreq),(1.0/nfreq))
    discarded_phones = []

    for binary_file in binary_data_files:

        print("Working on file %s" %binary_file)

        cell_ID = binary_file[10:-4]
        metadata_file = 'metadata_cell_'+cell_ID+'.pkl'

        metadata = pickle.load(open(metadata_file, "rb"))
        binfile = open(binary_file,'rb')
        df = pd.DataFrame(metadata)

        #Generate a new time vector on which to interpolate all values
        cumulative_sum = np.concatenate((np.array([0]),np.array(df['nbytes_offset'].cumsum().values)))
        offsets = df['nbytes_offset'].values

        dphones = 0

        for i in range(len(df)):

            mdata = df.iloc[i]
            #This is a pandas timestamp object
            start_time = mdata['tstart']
            end_time = mdata['tend']
            duration = (end_time-start_time).total_seconds()

            #only proceed if we've recorded more than 5 mins of data

            if duration > nsec:

                nb_offset = cumulative_sum[i]
                nb_read = offsets[i] 

                binfile.seek(0)
                binfile.seek(nb_offset)
                data=binfile.read(nb_read)
                nend = int(len(data)/4)
                data=np.asarray(struct.unpack(nend*'f',data[:]))
                data=np.reshape(data,(4,int(len(data)/4)))[:,0:nend]

                X = data[0]
                Y = data[1]
                Z = data[2]
                T = data[3]

                try:
                    resample(X, Y, Z, T, Tnew)
                except:
                    print("Error in resample")
                    dphones += 1 

            else:

                print("Record length is < 5 mins!")
                dphones += 1

        discarded_phones.append(dphones)


if __name__ == '__main__':

    t0 = time.time()
    main()
    t1 = time.time()
    print('Total time: %f' %(t1-t0))

