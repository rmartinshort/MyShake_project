#!/usr/bin/env python

#RMS 2018
#Adapted from QK. 
#Search MongoDB for triggers. We can use information from this to report the two-way trigger time
#This may be a more accurate way of determining trigger delay than using the times from the 
#sensordata database

from datetime import datetime
import pytz
import pandas as pd
import numpy as np
from datetime import timedelta, datetime
import pickle
import smtplib
from pymongo import MongoClient, ReadPreference 
import time
import os


def read_triggers_mongoDB(password_file, t0, t1, ids, tz = 'UTC', port = 27017):

    #Start and end times of the recording as strings
    t0 = t0.strftime("%Y-%m-%d %H:%M:%S") 
    t1 = t1.strftime("%Y-%m-%d %H:%M:%S")

    '''Read from the mongoDB between t0 and t1. All devices with id in ids'''

    def timestamps2str(timestamps, tz = 'UTC'):
        if tz == 'UTC':
            df_tmp = timestamps.apply(lambda ts: datetime.fromtimestamp(ts/1000., tz=pytz.utc))
        else:
            df_tmp = timestamps.apply(lambda ts: datetime.fromtimestamp(ts/1000., tz=pytz.timezone('US/Pacific')))    
        return df_tmp

    def str2timestamps(timeStr):

        try:
            return timeStr.apply(lambda s: (datetime.strptime(s.strip(), '%Y-%m-%d %H:%M:%S') - datetime(1970, 1, 1)).total_seconds())
        except:    
            return timeStr.apply(lambda s: (datetime.strptime(s.strip(), '%Y-%m-%d %H:%M:%S.%f') - datetime(1970, 1, 1)).total_seconds())


    identity = pickle.load(open(password_file, 'rb'))
    c = MongoClient('bsl-myshake01', port)
    c['eew'].authenticate(identity['mongoUser'], identity['mongoPasswd'])
    
    t_start = time.mktime(datetime.strptime(t0, "%Y-%m-%d %H:%M:%S").timetuple())*1000
    t_end = time.mktime(datetime.strptime(t1, "%Y-%m-%d %H:%M:%S").timetuple())*1000
    
    query_results = list(c['eew'].triggers.find({"$and" : [{"tt":{"$lt": t_end, "$gt": t_start}}, {"deviceid":{"$in":ids}}, {'tf':{"$eq":6}} ]}))
 
    df = pd.DataFrame(query_results)
    print(df)
    df['triggerTime'] = timestamps2str(df['tt'], tz = tz)
    df = df.set_index('triggerTime')

    return df
  

##### For testing

# t0 = '2018-09-01 00:00:00'
# t1 = '2018-09-02 00:00:00'

# ids = ["b8060d06-fdce-31e3-a497-4f949414c42e",
# "e3172a2f-9d0a-3bc1-83cd-9520b5c93372",
# "527184f7-dd7a-3251-a23f-7d9720e2d7b2",
# "dce6b393-cba8-3fc3-9630-91ed3ff2acc5"]

# df = read_triggers(t0, t1, ids, tz = 'UTC', port = 27017)
# print(df)
# df.to_csv('test1.csv')

#time_range = t0.replace(' ' , 'T').replace(':', '').replace('-', '') + '_' + t1.replace(' ' , 'T').replace(':', '').replace('-', '')
#pickle.dump(df, open('./' + time_range + '_triggers.pkl', 'wb'))

