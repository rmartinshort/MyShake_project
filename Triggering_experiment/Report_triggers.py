#!/usr/bin/env python

#AI 2018. Modified by RMS 2018
#Modified for clarity and to work on bsl-myshake-rnd01

#This should be run as a check to ensure that phones have been correclly been triggered. It looks at the waveform database
#and reports phones within the region of interest that have been activated in response to the triggering command

import numpy as np
import os
import json
import pickle
from datetime import timedelta, datetime
import matplotlib.pyplot as plt
import sys
import shutil
import ssl
import glob
import gzip
from pandas import HDFStore, read_hdf, Timedelta
import pytz

#Set path to Asaf's functions
sys.path.append('/work/myshake/rmartin/asaf_scripts/get_data/functions/')
#import from Asaf's functions
from functions_io import *



def buildSQL_sliceTime(t_b, t_e, lon_min, lon_max, lat_min, lat_max):

    '''Build SQL query to extract triggers that occured between t_b and t_e within some region'''
    
    t_1 = t_b.strftime("%Y-%m-%d %H:%M:%S") 
    t_2 = t_e.strftime("%Y-%m-%d %H:%M:%S") 
    sql = "select * from waveform where triggerFrom=6 and latitude < " + str(lat_max) + " and latitude > " + str(lat_min) + " and longitude < " + str(lon_max) + " and longitude > " + str(lon_min) + " and triggerTimer between UNIX_TIMESTAMP(\'" + t_1 + "\')*1000 and UNIX_TIMESTAMP(\'" + t_2 + "\')*1000;"

    return sql

def get_hb(HBfile):

    '''Get all the heartbeat data from provided file in the current directory. Generate a dataframe'''

    df_s = []

    tb=datetime.strptime(HBfile.split('_')[3]+' '+HBfile.split('_')[4].split('-')[0],"%Y-%m-%d %H:%M:%S")
    te=datetime.strptime(HBfile.split('_')[4][9:19]+' '+ HBfile.split('_')[5].split('U')[0],"%Y-%m-%d %H:%M:%S")
        
    df_s.append(pickle.load(open(HBfile,'rb')))
    
    df_s=pd.concat(df_s) 
    df_s['ts']=pd.to_datetime(df_s['ts'],unit='ms')

    return df_s, tb, te

def ssl_connect(t_b,t_e,lat_min,lat_max,lon_min,lon_max,password_file):

    '''
    Apply the SQL query
    '''
     
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    identity = pickle.load(open(password_file, "rb"))

    #Connect to the sensordata database
    # Remember the port number: 3307
    engine_wf = create_engine('mysql+mysqldb://' + identity['user2'] + ':' + identity['passwd2'] + '@bsl-myshake01.geo.berkeley.edu:3307/sensordata')
    conn_wf = engine_wf.connect()
    sql = buildSQL_sliceTime(t_b,t_e,lon_min,lon_max,lat_min,lat_max)
    df_from_sql = pd.read_sql_query(sql, conn_wf)
    df_from_sql['triggerTimer']=pd.to_datetime(df_from_sql['triggerTimer'],unit='ms')

    return df_from_sql


def main(region_code='LA'):

    '''
    Set up region information. Extract the phones that have been triggered and display their information
    '''

    #We want to convert everything to UTC
    global local_tz, utc_tz
    local_tz = pytz.timezone("US/Pacific")
    utc_tz = pytz.timezone("utc")
    password_file = 'ttt.identity.p'

    #regions = {'LA':[33,34.5,-118.7,-117.2]}
    regions = {'LA':[33,34.5,-118.7,-117.2],'SFbay':[36.921,38.494,-123.128,-121.387]}

    #Set up the coordinates of the region we wish to investiage
    region = regions[region_code]
    lat_min = region[0]
    lat_max = region[1]
    lon_min = region[2]
    lon_max = region[3]

    #The region we are trying to pull phones from
    print('Region: %s [%s/%s/%s/%s]' %(region_code,lat_min,lat_max,lon_min,lon_max))

    #This is the file from which to extract device information
    #HBfile = 'HB_reg_LA_2018-09-07_12:37:56-2018-09-07_14:37:56UTC.pkl'
    #Triggered phones
    #Triggered_phones = 'Phones_to_trigger_1536288122.csv'
    #Planned_triggered_phones = pd.read_csv(Triggered_phones)
    #List of the phones that we aimed to trigger
    #Planned_triggered_ids = Planned_triggered_phones['deviceId']


    #Get dataframe of heartbeat info between t0 and t1, for the region of interest
    #The HBfile is generated by 'Extract_phones_to_trigger.py'
    #Note that t1 is the time when the recording was initiated

    #heartbeat_df, t0, t1 = get_hb(HBfile)

    #print(t0,t1)

    #We're interested in all the triggers between the start of the recording and some
    #small time afterwards.

    #t1 = datetime(2018,9,6)
    #t2 = datetime(2018,9,8)

    #t2 = t1 + timedelta(seconds=600)

    #print(t1,t2)

    t1 = datetime(2018,7,1)
    t2 = datetime(2018,9,10)

    #tb_all=local_tz.localize(t0,is_dst=None).astimezone(pytz.utc)
    #te_all=local_tz.localize(t2,is_dst=None).astimezone(pytz.utc)

    #Get dataframe of triggered phones
    all_triggered_phones = ssl_connect(t1,t2,lat_min,lat_max,lon_min,lon_max,password_file)

    #Get a list of the devices that triggered, ready for querying the mongodb
    #devices_list = list(all_triggered_phones['deviceId'].values)
    #mongodb_triggers = read_triggers_mongoDB(password_file,t1,t2,devices_list)

    #print(mongodb_triggers)

    #Select only the phones that we planned to trigger
    #all_triggered_phones = all_triggered_phones[all_triggered_phones['deviceId'].isin(Planned_triggered_ids)]

    #Write triggered phones to a file 
    print(all_triggered_phones)

    print('All triggered phones in Region: %s [%s/%s/%s/%s]' %(region_code,lat_min,lat_max,lon_min,lon_max))
    print('Between t1 = %s and t2 = %s' %(t1,t2))

    all_triggered_phones.to_csv("Devices_triggered_t1_t2.csv")


if __name__ == '__main__':

    main(region_code='LA')
