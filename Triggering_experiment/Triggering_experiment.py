#!/usr/bin/env python
 
#RMS 2018

#We want to invesigate the ability of triggering commands to actually be sent to phones in the network and the data that they send 
#back. To do this, we will send the trigger command to phones in phases over some length of time and moniter how many 
#end up actually being triggered.

import numpy as np
import pickle
from datetime import timedelta, datetime
import sys
import ssl
from pandas import HDFStore, read_hdf, Timedelta
import pytz
import time

#Set path to Asaf's functions
sys.path.append('/work/myshake/rmartin/asaf_scripts/get_data/functions/')
#import from Asaf's functions
from functions_io import *

    
def buildSQL_sliceTime_hbHistory(t0, t1,lat_min,lat_max,lon_min,lon_max):

    ''' 
    This builds an SQL query for the database of phone information to extract just the phones within the selected time 
    and space frame
    '''

    sql = 'select ts,deviceId, latitude, longitude, hbSource,appVersion,profileName,profileVersion from hbHistory where ts<UNIX_TIMESTAMP(\'' + t1 + '\')*1000 and ts>UNIX_TIMESTAMP(\'' + t0 + '\')*1000 \
and latitude> '+ str(lat_min) +' and  latitude< '+ str(lat_max) +' and longitude > ' + str(lon_min) + ' and longitude < ' + str(lon_max) + ';'
    return sql

def buildSQL_sliceTime_waveform(t_b, t_e, lon_min, lon_max, lat_min, lat_max):

    '''Build SQL query to extract triggers that occured between t_b and t_e within some region'''
    
    t_1 = t_b.strftime("%Y-%m-%d %H:%M:%S") 
    t_2 = t_e.strftime("%Y-%m-%d %H:%M:%S") 
    print(t_1,t_2)
    sql = "select deviceId, triggerTimer, latitude, longitude, accuracy, altitude, triggerFrom, createdOn from waveform where triggerFrom=6 and latitude < " + str(lat_max) + " and latitude > " + str(lat_min) + " and longitude < " + str(lon_max) + " and longitude > " + str(lon_min) + " and triggerTimer between UNIX_TIMESTAMP(\'" + t_1 + "\')*1000 and UNIX_TIMESTAMP(\'" + t_2 + "\')*1000;"

    return sql

def ssl_connect_hbHistory(password_file,region_code,t0,t1,lat_min,lat_max,lon_min,lon_max):

    '''
    Apply sql query to hbHistory database
    '''

    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    identity = pickle.load(open(password_file, "rb"))

    # Connect to the myshake database
    # Remember that we need to use port 3306
    #engine_wf = create_engine('mysql+pymysql://' + 'myshake_research' + ':' + identity['passwd1'] + '@bsl-myshake01/myshake?unix_socket=' + identity['socket1'])
    engine_wf = create_engine('mysql+mysqldb://' + identity['user1'] + ':' + identity['passwd1'] + '@bsl-myshake01.geo.berkeley.edu:3306/myshake')
    conn_wf = engine_wf.raw_connection()
    sql = buildSQL_sliceTime_hbHistory(t0, t1,lat_min,lat_max,lon_min,lon_max)
    df_tmp = pd.read_sql_query(sql, conn_wf)
    pickle.dump(df_tmp,open('HB_reg_%s_%s_%s-%s_%sUTC.pkl'%(region_code,t0.split()[0],t0.split()[1],t1.split()[0],t1.split()[1]),'wb'))
    return df_tmp

def ssl_connect_waveform(t_b,t_e,lat_min,lat_max,lon_min,lon_max,password_file):

    '''
    Apply the SQL query to waveform database 
    '''
     
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    identity = pickle.load(open(password_file, "rb"))

    #Connect to the sensordata database
    # Remember the port number: 3307
    engine_wf = create_engine('mysql+mysqldb://' + identity['user2'] + ':' + identity['passwd2'] + '@bsl-myshake01.geo.berkeley.edu:3307/sensordata')
    conn_wf = engine_wf.connect()
    sql = buildSQL_sliceTime_waveform(t_b,t_e,lon_min,lon_max,lat_min,lat_max)
    df_from_sql = pd.read_sql_query(sql, conn_wf)
    df_from_sql['triggerTimer']=pd.to_datetime(df_from_sql['triggerTimer'],unit='ms')
    return df_from_sql

def Trigger_phones(password_file,trigger_script,region_code,t0_dt,t1_dt,lat_min,lat_max,lon_min,lon_max):

    '''
    Write a file containing the phones we want to trigger, then go ahead and actually trigger them
    '''

    outfilename = 'Phones_to_trigger.txt'

    #Start and end times of the recording
    t0 = t0_dt.strftime("%Y-%m-%d %H:%M:%S") 
    t1 = t1_dt.strftime("%Y-%m-%d %H:%M:%S") 

    print("Querying database for phone heartbeats in region of interest")

    df_tmp=ssl_connect_hbHistory(password_file,region_code,t0,t1,lat_min,lat_max,lon_min,lon_max)

    print("Done query")

    #Convert hbtime to time units
    df_tmp['hbtime']=pd.to_datetime(df_tmp['ts'],unit='ms')
    tmp=df_tmp.drop_duplicates(['deviceId'])
    #ensure that users have the latest app
    tmp = tmp[tmp['appVersion'].values.astype(np.int32)>200]

    #Write the test_dev.txt file
    fp=open(outfilename,'w')
    fp.write("body={\"deviceList\": [")

    for i in range(len(tmp['deviceId'])):
        fp.write("\"%s\""%(tmp['deviceId'].values[i]))

        if i<len(tmp['deviceId'])-1:
            fp.write(",")

    fp.write('],"path": "AsafTest"}')
    fp.close()

    #Send trigger message
    print ("Sending trigger message to phones")
    os.system('./%s' %trigger_script)
    trigger_send_time = datetime.now()

    #Write the triggered device dataframe to a .csv file
    attempted_triggers = "Phones_to_trigger_%s.csv" %trigger_send_time.strftime('%s')
    tmp.to_csv(attempted_triggers)

    return attempted_triggers, trigger_send_time


def main(region_code='LA'):

    '''
    Set the region information and generate list of phones that had sent HB data in the last nhours hours. This can then 
    be used to extract waveform information from these phone (i.e. turn them on)
    '''

    trigger_script = 'rec_start_rms.sh'

    regions = {'LA':[33,34.5,-118.7,-117.2],'SFbay':[36.921,38.494,-123.128,-121.387]}

    #Set up the coordinates of the region we wish to investiage
    region = regions[region_code]
    lat_min = region[0]
    lat_max = region[1]
    lon_min = region[2]
    lon_max = region[3]

    #The region we are trying to pull phones from
    print('Region: %s [%s/%s/%s/%s]' %(region_code,lat_min,lat_max,lon_min,lon_max))

    #Experiment set up
    local_tz = pytz.timezone("US/Pacific")
    t0_start = datetime.today() #start time of the experiment 
    runtime = 48*3600 #length of time to run the experiment
    trigger_interval = 40*60 #how often we will send command to trigger new phones
    nhours=1 #This is the number of hours before the search was launched in which to search for phones that sent HB data
    endtime = t0_start + timedelta(seconds=runtime)
    currenttime = t0_start

    print("Experiment: starttime: %s  endtime: %s  trigger_interval: %s" %(t0_start,endtime,trigger_interval))

    #path to ttt.identity.p. Put this file somewhere on your home dir. at myshake01
    password_file="/work/myshake/rmartin/asaf_scripts/ttt.identity.p"

    while currenttime <= endtime:

        t0_dt = currenttime + timedelta(hours=nhours*-1) #start time of search
        t1_dt = currenttime #end time of search (now)

        #t0_local = t0_dt.replace(tzinfo=pytz.utc).astimezone(local_tz)
        #t1_local = t1_dt.replace(tzinfo=pytz.utc).astimezone(local_tz)
        t0_local = t0_dt
        t1_local = t1_dt
        
        print ('Start time: %s. End time: %s' %(t0_local,t1_local))

        #Start recording on these phones
        attempted_triggers, trigger_send_time = Trigger_phones(password_file,trigger_script,region_code, t0_local, t1_local, lat_min, lat_max, lon_min, lon_max)

        print('Trigger send time: %s' %trigger_send_time)

        Planned_triggered_phones = pd.read_csv(attempted_triggers)

        #List of the phones that we aimed to trigger
        Planned_triggered_ids = Planned_triggered_phones['deviceId']

        #It seems to take quite a while for triggered phones to show up in the database
        #We need to wait some time before we can proceed to see what has been triggered
        time.sleep(20*60)

        #Moniter the recording on thse phones between the send time and 10 minutes after the send time
        t2_local = trigger_send_time + timedelta(seconds=600)
        #Get dataframe of triggered phones during the time frame and area that we are interested in
        all_triggered_phones = ssl_connect_waveform(trigger_send_time,t2_local,lat_min,lat_max,lon_min,lon_max,password_file)

        #Select only the phones that we planned to trigger
        all_triggered_phones = all_triggered_phones[all_triggered_phones['deviceId'].isin(Planned_triggered_ids)]
        print(all_triggered_phones)

        all_triggered_phones.to_csv("Devices_triggered_%s_%s.csv" %(trigger_send_time.strftime('%s'),t2_local.strftime('%s')))

        #Wait until its time to trigger again
        time.sleep(trigger_interval)

        #Update the time counter
        currenttime = datetime.now()

if __name__ == '__main__':
    main()

        
         
