#!/usr/bin/env python

#AI 2018. Modified by RMS 2018
#Write all triggers within some time period to a hdf file

import numpy as np
import os
import json
import pickle
from datetime import timedelta, datetime
from obspy.geodetics.base import gps2dist_azimuth
import sys
import shutil

#This is a way of accessing code stored in the functions directory
sys.path.append('./functions')
from functions_io import *
import ssl
import gzip

hd5='myshakeMeta.h5'
FID='ttt.identity.p'

#Writing all the triggered phone data to a hdf5 file

try:
    store = pd.HDFStore(hd5,mode='r+')
except:
    store = pd.HDFStore(hd5,mode='w')
#t0=datetime(2016,2,20,0,0,0)
t0=datetime(2017,9,20,8,53,14)
t1=datetime.now()
ndays=2
while t0<t1:
    
    print (t0)
    t0_tmp = t0.strftime("%Y-%m-%d %H:%M:%S")
    t1_tmp = (t0+timedelta(days=ndays)).strftime("%Y-%m-%d %H:%M:%S")


    identity = pickle.load(open(FID, "rb"))

    #Generate connection to the database
    engine = create_engine('mysql+pymysql://' + identity['user1'] + ':' + identity['passwd1'] + '@localhost/myshake?unix_socket=' + identity['socket1'])
    conn = engine.connect()

    #From functions_io.py. Make a dataframe of the heartbeat information
    df_heartbeat=read_heartbeat(t0_tmp, t1_tmp, conn)
    
    del df_heartbeat['appVersion'],df_heartbeat['annVersion']
    store.append('Heartbeat', df_heartbeat,data_columns=True)

    #From functions_io.py
    df_trig=read_triggers(t0_tmp, t1_tmp, FID=FID)

    longitude=(pd.DataFrame(df_trig['l'].apply(pd.Series)['coordinates'].apply(pd.Series))).values[:,0].copy()
    latitude=(pd.DataFrame(df_trig['l'].apply(pd.Series)['coordinates'].apply(pd.Series))).values[:,1].copy()
    df_trig['_id']=df_trig['_id'].astype(str)
    df_trig['deviceid']=df_trig['deviceid'].astype(str)
    df_trig['longitude']=longitude
    df_trig['latitude']=latitude
    del df_trig['l'],df_trig['schema_version']

    store.append('Trigger', df_trig,data_columns=True)
    t0+=timedelta(days=ndays)
    break

store.close()

