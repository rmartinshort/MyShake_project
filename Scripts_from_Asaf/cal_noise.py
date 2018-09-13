#!/usr/bin/env python


from pandas import HDFStore, read_hdf, Timedelta
import pandas as pd
import matplotlib.pyplot as plt
import collections
import numpy as np
import math
import time
from datetime import datetime,timedelta
from geopy.distance import vincenty
from pygc import great_circle
import pytz
#import pickle
import _pickle as pickle
from itertools import count
from numba import jit
import multiprocessing
from multiprocessing import Manager
from multiprocessing.sharedctypes import Value, Array
from ctypes import c_double,c_int
import gzip
import os
from obspy import Trace,Stream
#from obspy.core.inventory import Inventory, Network, Station, Channel, Site
#import obspy
import h5py
import matplotlib.path as mplPath
#import pyasdf
from itertools import *
import sys

def read_cat(region):

    fp=open('../aux/Cat_%s.pyc'%(region),'rb')
    Cat=pickle.load(fp)
    u, indices = np.unique(np.asarray([t for t in Cat['eve_id']]), return_index=True) 
    key_del=-1
    print (Cat.keys())
    for key in Cat.keys():
        if len(Cat[key])==0:
            key_del=key
            continue
        Cat[key]=[Cat[key][i] for i in indices]
    fp.close()
    if key_del!=-1:
        del Cat[key_del]
    print (Cat.keys())
    return Cat

def get_dist(sta_lon,sta_lat,eve_lon,eve_lat):

    sta_lon=np.asarray(sta_lon)
    sta_lat=np.asarray(sta_lat)
    dist=[]

    for slat,slon in zip(sta_lat,sta_lon):
        dist.append(vincenty((eve_lat,eve_lon),(slat,slon)).kilometers)

    return np.asarray(dist)



def read_cat_associate_trig(df_trig,region):
    Cat=read_cat(region)
    
    neve_in_cat=len(Cat['origin_time'])
    Output={}
    Output['eve']=[]
    Output['epi_dist']=[]
    Output['potential_trig_15km_dist']=[]
    Output['df_trig_potent']=[]
    mo0=Cat['origin_time'][0].month

    for ieve in range(neve_in_cat):
        t_eq=Cat['origin_time'][ieve]
        mask=((df_trig.index >= t_eq) & (df_trig.index <= t_eq+timedelta(seconds=120)))
        df_trig_potential=df_trig[mask]
        if len(df_trig_potential)==0:
            continue
        
        dist=get_dist(df_trig_potential['longitude'],df_trig_potential['latitude'],Cat['lon'][ieve],Cat['lat'][ieve])
        indx_sta_in_rad=np.where(dist<=15)[0]
        
        print (t_eq,len(df_trig_potential),Cat['mag'][ieve],np.min(dist),np.max(dist))
        Output['eve'].append({key : Cat[key][ieve] for key in Cat.keys()})
        Output['epi_dist'].append(dist)
        Output['potential_trig_15km_dist'].append(indx_sta_in_rad)
        Output['df_trig_potent'].append(df_trig_potential)

        if t_eq.month!=mo0:
            print (t_eq.month)
            mo0=t_eq.month

    return Output

#function to get ntrigs from object output from read_cat_associate_trig
#region = Ncal,Scal
#epi_dist_max = search radius
def slice_df_by_epiDist(Obj,epi_dist_max):
    
    Output={}
    Output['eve']=[]
    Output['epi_dist']=[]
    Output['potential_trig_dist']=[]
    Output['df_trig_potent']=[]

    for ieve,eve in enumerate(Obj['eve']):

        df_trig_potential=Obj['df_trig_potent'][ieve]
        eve_lon=eve['lon']
        eve_lat=eve['lat']
        dist=get_dist(df_trig_potential['longitude'],df_trig_potential['latitude'],eve_lon,eve_lat)
        indx_sta_in_rad=np.where(dist<=epi_dist_max)[0]
        if len(indx_sta_in_rad)==0:
            continue

        Output['eve'].append(eve)
        Output['epi_dist'].append(dist)
        Output['potential_trig_dist'].append(indx_sta_in_rad)
        Output['df_trig_potent'].append(df_trig_potential)

    mag,ntrig = plot_mag_detect(Output)

    return (mag,ntrig)

def get_events_from_client(reg):
    
    from obspy.clients.fdsn import Client
    from obspy.core.utcdatetime import UTCDateTime

    win=get_spat_win(reg)

    
    if reg is 'Ncal':
        client = Client("NCEDC")
    else:
        #client = Client("SCEDC")
        tmp={}
        tmp['eve']=[]

        fp=open('../aux/cat_Scal_201602-201709.txt','r')
        for line in fp:
            eve={}
            eve['eventid']=line.split()[0]
            eve['event_type']=line.split()[1]
            eve['lat']=float(line.split()[4])
            eve['lon']=float(line.split()[5])
            eve['depth']=float(line.split()[6])
            eve['mag']=float(line.split()[7])
            T=datetime.strptime(line.split()[3].split('.')[0]+'.'+line.split()[3].split('.')[1]+'000','%Y/%m/%d,%H:%M:%S.%f')
            eve['origin_time']=T
            tmp['eve'].append(eve)


        fp.close()
        pickle.dump(tmp,open('Cat_%s.pyc'%(reg),'wb'))

        return


    starttime = UTCDateTime("2016-02-01")
    endtime = UTCDateTime(datetime.now())
    #endtime=UTCDateTime("2016-02-07")

    print (win)
    
    cat = client.get_events(starttime=starttime, endtime=endtime,minlatitude=win['lat'][0],maxlatitude=win['lat'][1],\
                            minlongitude=win['lon'][0],maxlongitude=win['lon'][1])
    
    tmp={}
    tmp['eve']=[]
    for e in cat:
        eve={}
        if (len(e['magnitudes'])==0) or e['magnitudes'][0]['mag']<1.5:
            continue
        try:
            eve['lat']=e['origins'][0]['latitude']
            eve['lon']=e['origins'][0]['longitude']
            eve['depth']=e['origins'][0]['depth']*1e-3
            eve['mag']=e['magnitudes'][0]['mag']
            if eve['lat']>38.68 and eve['mag']<2.5:
                continue
            
            
            eve['origin_time']=e['origins'][0]['time']
            eve['eventid']=e['extra']['eventid']['value']
            eve['event_type']=e['event_type']
        except:
            continue
        tmp['eve'].append(eve)

    pickle.dump(tmp,open('Cat_%s.pyc'%(reg),'wb'))



#find potential traces around the origin time from dictionary Obj 
def get_rec_time_space_window(Obj,spatial_win):
    
    sec_before=timedelta(seconds=240)
    sec_after=timedelta(seconds=240)
    lat_min=spatial_win['lat'][0]
    lat_max=spatial_win['lat'][1]
    lon_min=spatial_win['lon'][0]
    lon_max=spatial_win['lon'][1]
    hd5="/home/ainbal/Data/myshakeMeta_before_20161205.h5"
    dir_data="/data/sensordata/output"
    Output={}
    Output['eve_files']={}
    Output['eve_files']['eve']=[]
    Output['eve_files']['files']=[]
    slat=37.8762 
    slon=-122.2356

    for ieve,eve in enumerate(Obj['eve']):
        t_eq=eve['origin_time']

        #NOTE
        dist=vincenty((eve['lat'],eve['lon']),(slat,slon)).kilometers
        if dist>80:
            continue
       
        print (ieve,eve)
        time_win=[t_eq-sec_before,t_eq+sec_after]
        qbeg='%d%02d%02d %02d:%02d:00'%(time_win[0].year,time_win[0].month,time_win[0].day,time_win[0].hour,time_win[0].minute)
        qend='%d%02d%02d %02d:%02d:00'%(time_win[1].year,time_win[1].month,time_win[1].day,time_win[1].hour,time_win[1].minute)
        query = 'index>\"%s\" & index<\"%s\"'%(qbeg,qend)
        
        hdf_trig = read_hdf(hd5, key = 'Trigger', where = query)
        
        hdf_trig = hdf_trig[['deviceid','tt', 'latitude', 'longitude', 'tf','ti']]
        hdf_trig = hdf_trig[(hdf_trig['latitude']>lat_min) & (hdf_trig['latitude']<lat_max) & (hdf_trig['longitude']>lon_min) &
                           (hdf_trig['longitude']<lon_max)]
        fnames=[]
        dev_list=[]
        for trig_time,dev,tt,ti in zip(hdf_trig.index,hdf_trig['deviceid'],hdf_trig['tt'],hdf_trig['ti']):
            if (len(hdf_trig['deviceid'])<10):
                break
            path=build_path(dir_data,dev,tt)
            if path==-1:
                continue
                print ('%s %s %lf'%(dev,tt,ti))
            fnames.append(path)
            dev_list.append(dev)
        
        if (len(dev_list)==0):
            continue

        u, indices = np.unique(np.asarray(dev_list), return_index=True)

        Output['eve_files']['eve'].append(eve)
        Output['eve_files']['files'].append([fnames[i] for i in indices])


    #sorting by number of potential traces
    n=[len(f) for f in Output['eve_files']['files']]
    indices=np.argsort(n)[::-1]
    print (indices)
    Output['eve_files']['eve']=[Output['eve_files']['eve'][i] for i in indices]
    Output['eve_files']['files']=[Output['eve_files']['files'][i] for i in indices]

    return Output

#buile path to data files
def build_path(dir_data,dev,tt):
    
    T=float('%s.%s'%(str(tt)[0:10],str(tt)[10:]))
    t0=datetime.fromtimestamp(T)
    year=t0.year
    month=t0.month
    day=t0.day
    hour=t0.hour

    path='%s/%d/%02d/%02d/%02d:00:00/%s_%s.pkl.gz'%(dir_data,year,month,day,hour,dev,tt)
    if os.path.isfile(path) is False:
        print (path)
        return -1
    return path

#read data files
def read_file_from_eve(Obj,indx):

    
    #NOTE for ncal events:
    #indx=6
    
    Output={}
    Output['trace']=[]
    event=Obj['eve_files']['eve'][indx]
    
    Output['event']=event
    print ('read traces for ',event)
    #Output['trace']['time']=[]
    #Output['trace']['x']=[]
    #Output['trace']['y']=[]
    #Output['trace']['z']=[]

    for fname in Obj['eve_files']['files'][indx]:
        try:
            inF = gzip.open(fname, "rb")
        except:
            print (fname)
            continue
        obj=pickle.load(inF,encoding='latin1')
        #devId=obj
        #tsamp=(pd.to_datetime(obj['data']['ts'],unit='ms'))
        Output['trace'].append(obj)
    
    return Output


#a function that groups tt of indevidual phones
#will output the number of detections for eahc phone, which will
#be used as cutoff criteria
def get_trig_dev(spatial_win):
    lat_min=spatial_win['lat'][0]
    lat_max=spatial_win['lat'][1]
    lon_min=spatial_win['lon'][0]
    lon_max=spatial_win['lon'][1]
    hd5="/home/ainbal/Data/myshakeMeta_before_20161205.h5"
    dir_data="/data/sensordata/output"
    time_win=[datetime(2016,2,10),datetime(2017,1,1)]
    ndays=(time_win[1]-time_win[0]).days
    Output={}
        
    qbeg='%d%02d%02d %02d:%02d:00'%(time_win[0].year,time_win[0].month,time_win[0].day,time_win[0].hour,time_win[0].minute)
    qend='%d%02d%02d %02d:%02d:00'%(time_win[1].year,time_win[1].month,time_win[1].day,time_win[1].hour,time_win[1].minute)
    query = 'index>\"%s\" & index<\"%s\"'%(qbeg,qend)
    
    hdf_trig = read_hdf(hd5, key = 'Trigger', where = query)
    hdf_trig = hdf_trig[['deviceid','tt', 'latitude', 'longitude', 'tf','ti']]
    hdf_trig = hdf_trig[(hdf_trig['latitude']>lat_min) & (hdf_trig['latitude']<lat_max) & (hdf_trig['longitude']>lon_min) &\
                        (hdf_trig['longitude']<lon_max)]

    hdf_trig['one']=np.ones(len(hdf_trig))
    
    g = pd.DataFrame({'count' : hdf_trig.groupby('deviceid').size()}).reset_index()
    return (g,ndays)
    #print (g['deviceid'].iloc)



def read_data_from_Traces(traces):

    from obspy import Trace,Stream
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.core.trace import Stats
    from scipy.interpolate import interp1d
   
    out=[]
    eve=traces['event']

    for comp in ['x','y','z']:
        arr=[]
        
        for t in traces['trace']:
            
            ttrig=pd.to_datetime(t['triggerTimer'],unit='ms')
            tsamp=(pd.to_datetime(t['data']['ts'],unit='ms'))
            torigin=eve['origin_time']

            if tsamp[len(tsamp)-1] < torigin:
                continue
            
            T=Trace()

            stats=Stats()
            stats.network='MyShake'
            stats.sampling_rate=25
            stlon=t['location']['longitude']
            stlat=t['location']['latitude']
            evelon=eve['lon']
            evelat=eve['lat']
            epidist=vincenty((evelat,evelon),(stlat,stlon)).kilometers
            if epidist>10:
                continue
            mag=eve['mag']

            stats.location=t['location']
            stats.channel=comp
            stats.station=t['deviceId']
            stats.evelat=evelat 
            stats.evelon=evelon
            stats.distance=epidist*1e3
            stats.torigin=torigin
            stats.mag=mag
            stats.ttrig=ttrig
           
            Dt=UTCDateTime(torigin)

            nsamp=len(t['data']['ts'])
            tbeg=t['data']['ts'][0]
            tend=t['data']['ts'][nsamp-1]

            new_time=np.arange(tbeg,tend,20)
            gaps=find_data_gaps(t['data']['ts'],nsamp)
            f=interp1d(t['data']['ts'],t['data'][comp])

            if len(gaps)==0:
                acc=f(new_time)

            else:
                tmp=np.zeros(len(np.arange(tbeg,tend,20)))
                i=0
                sign=0
                for g in gaps:
                    if (g[2] > torigin and g[3] < torigin+timedelta(seconds=50)) or\
                       (g[2] < torigin and g[3]>torigin ) :
                        print (g[2],g[3],torigin,torigin+timedelta(seconds=50))
                        sign=1
                        break

                    iend=g[0]
                    #print (i,iend)
                    tmp[i:iend]=f(new_time[i:iend])
                    i=g[1]
                
                #if len(gaps)>2:
                #    for g in gaps:
                #        print (tmp[g[0]:g[1]])
                #    return
                if sign==1:
                    continue
                acc=tmp
            
            stats.starttime=UTCDateTime(tsamp[0])
            T.stats=stats

            T.data=acc
            T=T.slice(Dt-40,Dt+40)
            if len(T.data)==0:
                continue
            if epidist>20:
                continue
            T=T.detrend('demean')
            T=T.filter('bandpass',freqmin=2,freqmax=4)
            
            #print (tsamp[len(tsamp)-1],stats.endtime,stats.npts,nsamp,len(new_time),len(indx[0]))
            
            arr.append(T)

        stream=Stream(traces=arr) 
        out.append(stream)
    
    for j in range(len(out)):
        dist=[i.stats['distance'] for i in out[j]]
        indx=np.argsort(dist)
        out[j]=Stream([out[j][i] for i in indx])
        #out[j] = sorted(out[j], key=lambda out[j]: out[j].stats['distance'])

    return (out[0],out[1],out[2])

#nsta is the smartphone index
#st contains all smartphone traces 
def get_NCEDC_data(st,nsta):

    from obspy.clients.fdsn import Client
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.geodetics.base import gps2dist_azimuth
    from obspy import Trace,Stream

    t0=st[nsta].stats.starttime
    t1=st[nsta].stats.endtime
    down_start=t0-(t0.microsecond)/1e6
    down_end=t1-(t1.microsecond)/1e6+1
    originTime=UTCDateTime(st[nsta].stats.torigin)

    tmp=pickle.load(open('../aux/sta_list_ncedc_2000-2017.pyc','rb'))
    sta=[]
    for s in tmp:
        starttime=np.min(np.asarray([t for t in s['starttime']]))
        endtime=np.max(np.asarray([t for t in s['endtime']]))
        if UTCDateTime(endtime)<t1 or UTCDateTime(starttime)>t0:
            continue
        #epicentral distance
        #s['epidist']=vincenty((st[nsta].stats['evelat'],st[nsta].stats['evelon']),(s['latitude'],s['longitude'])).kilometers
        tmp2=gps2dist_azimuth(st[nsta].stats['evelat'],st[nsta].stats['evelon'],s['latitude'],s['longitude'])
        s['epidist']=tmp2[0]/1e3
        s['back_azimuth']=tmp2[2]
        #smartphone distance
        s['smdist']=vincenty((st[nsta].stats['location']['latitude'],st[nsta].stats['location']['longitude']),(s['latitude'],s['longitude'])).kilometers
        sta.append(s)

    del tmp
    

    #sort by epicentral dist.
    sta = sorted(sta, key=lambda sta: sta['epidist'])
    
    print (sta[0])
    client = Client("NCEDC")
    O_ep=[]
    for s in sta:
        
        for ich,chan in enumerate(s['channel']):
            
            if chan[1] not in 'HLGMNP' or chan[0] in 'LVURPE' :
                continue
            
            try:
                u=s['unit'][ich]
            except:
                break

            if UTCDateTime(s['endtime'][ich])<t1 or UTCDateTime(s['starttime'][ich])>t0:
                continue
            
            #print (s['net'], s['name'], s['loc'] , chan, down_start, down_end)
            try:
                struct = client.get_waveforms(network=s['net'], station=s['name'], location='*', channel=chan, starttime=down_start, endtime=down_end)
            
            except:
                #print ('no data for ',s['net'], s['name'],chan,)
                continue
                
            try: 
                print (s['net'], s['name'], "00", chan, down_start, down_end,s['unit'][ich])
            except:
                print (chan,len(O_ep),down_start, down_end,ich)
                return
            struct[0].stats['epidist']=s['epidist']
            struct[0].stats['smdist']=s['smdist']
            struct[0].stats['back_azimuth']=s['back_azimuth']
            struct[0].stats['longitude']=s['longitude']
            struct[0].stats['latitude']=s['latitude']
            struct[0].stats['torigin']=originTime
            
            if 'VEL' in s['unit']:
                struct=struct.differentiate()
             
            O_ep.append(struct[0])
            if len(O_ep)==3:
                break
        
        if len(O_ep)>0:
            break
        
        else:
            continue
    
    O_ep=Stream(O_ep)
    O_ep=O_ep.detrend('demean')
    O_ep=O_ep.filter('bandpass',freqmin=2,freqmax=4)
    #sort by smart dist.
    
    O_sm=[]
    sta = sorted(sta, key=lambda sta: sta['smdist'])
    
    for s in sta:

        for ich,chan in enumerate(s['channel']):
            
            if chan[1] not in 'HLGMNP' or chan[0] in 'SLVURPE' :
                continue
            
            if UTCDateTime(s['endtime'][ich])<t1 or UTCDateTime(s['starttime'][ich])>t0:
                continue
            
            #print (s['net'], s['name'], s['loc'] , chan, down_start, down_end)
            try:
                struct = client.get_waveforms(network=s['net'], station=s['name'], location='*', channel=chan, starttime=down_start, endtime=down_end)
            
            except:
                #print ('no data for ',s['net'], s['name'],chan,)
                continue
                
            
            print (s['net'], s['name'], "00", chan, down_start, down_end,s['unit'][ich])
            struct[0].stats['epidist']=s['epidist']
            struct[0].stats['smdist']=s['smdist']
            struct[0].stats['back_azimuth']=s['back_azimuth']
            struct[0].stats['longitude']=s['longitude']
            struct[0].stats['latitude']=s['latitude']
            struct[0].stats['torigin']=originTime
            
            if 'VEL' in s['unit']:
                struct=struct.differentiate()

            O_sm.append(struct[0])
            if len(O_sm)==3:
                break
        
        if len(O_sm)>0:
            break
        
        else:
            continue
    
    O_sm=Stream(O_sm)
    O_sm=O_sm.detrend('demean')
    O_sm=O_sm.filter('bandpass',freqmin=2,freqmax=4)
    print (O_ep[0].stats,O_sm[0].stats)

    return (O_ep,O_sm)



def find_data_gaps(T,nsamp):
    
    if nsamp%2!=0:
        T=np.append(T,T[nsamp-1]+20)
        nsamp+=1
    
    
    t0=T[0:nsamp:2]
    t1=T[1:nsamp:2]
    
    #print (len(T),len(t0),len(t1))
    
    DT=t1-t0
    indx=np.where(DT/1e3>1.0)
    gaps=[]
    for i in indx[0]:
        gaps.append((i*2,i*2+int(np.floor(DT[i]/20)),pd.to_datetime(t0[i],unit='ms'),pd.to_datetime(t1[i],unit='ms')))
        #print (nsamp,i*2,i*2+int(np.floor(DT[i]/20)))
    #print ()
    return gaps



def get_spat_win(region):

    #NOTE THAT I USED Scal HERE
    if region is 'LA':
        lat_min=33.35
        lat_max=34.1
        lon_min=-118.7
        lon_max=-117.4

    if region is 'Ncal':
        lat_min=35.455
        lat_max=40.3
        lon_max=-119
        lon_min=-127.2

    if region is 'SFBay':
        lat_min=36.921
        lat_max=38.494
        lon_min=-123.128
        lon_max=-121.387

    if region is 'Scal':
        lat_min=30.27
        lat_max=35.455
        lon_min=-120.5
        lon_max=-112

        
    
    spatial_win={}
    spatial_win['lat']=[lat_min,lat_max]
    spatial_win['lon']=[lon_min,lon_max]
    spatial_win['reg']=region
    
    
    return spatial_win

def read_dict_trig_epi_dist(region):
    fp=open('../aux/%s_obj_sta_trig_eve.pyc'%(region),'rb')
    Obj=pickle.load(fp)
    fp.close()
    return Obj


def plot_mag_detect(Obj):

    mag=[]
    ntrig=[]
    
    for i,indx_sta in enumerate(Obj['potential_trig_dist']):

        if len(indx_sta)==0:
            continue

        mag.append((Obj['eve'][i]['mag']))
        ntrig.append(len(indx_sta))

    return (mag,ntrig)

def dmean_filter(out,indx):
    
    from obspy import Trace,Stream
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.core.trace import Stats
    from scipy.interpolate import interp1d

    o0=out[0]
    o1=out[1]
    o2=out[2]

    s0=Stream([o0[i].detrend(type='demean') for i in range (indx,indx+1)])
    s1=Stream([o1[i].detrend(type='demean') for i in range (indx,indx+1)])
    s2=Stream([o2[i].detrend(type='demean') for i in range (indx,indx+1)])

    s0=s0.filter('bandpass',freqmin=2,freqmax=4)
    s1=s1.filter('bandpass',freqmin=2,freqmax=4)
    s2=s2.filter('bandpass',freqmin=2,freqmax=4)

    t=Stream([s0[0],s1[0],s2[0]])
    print (s0[0].stats.epidist,s0[0].stats.ttrig,s0[0].stats.torigin)
    t.plot(starttime=UTCDateTime(s0[0].stats.torigin),endtime=UTCDateTime(s0[0].stats.torigin)+20)


def Read_cat_2014():

    from obspy import read
    from obspy import Trace,Stream
    from obspy.core.utcdatetime import UTCDateTime
    import os

    fp=open('../aux/Cat_BayArea_2014.pyc','rb')
    Cat=pickle.load(fp)
    slat=37.8762 
    slon=-122.2356
    sta='BKS'
    ieve=0
    SM_chan=['DPX','DPY','DPZ']
    BK_chan=['BHE','BHN','BHZ']
    fmin=1
    fmax=2

    SM_dir='../Data/SmartphoneConRec2014/'
    BK_dir='../Data/Local_stations/'
    for eve_lon,eve_lat,eve_mag in zip(Cat['lon'],Cat['lat'],Cat['mag']):
        dist=vincenty((eve_lat,eve_lon),(slat,slon)).kilometers
        if eve_mag>2.2 and dist<10:
            obs_time=UTCDateTime(Cat['origin_time'][ieve])
            arr=[]
            for file in os.listdir(SM_dir):
                st=read('%s/%s'%(SM_dir,file),headonly=True)
                if st[0].stats.channel=='DPX' and st[0].stats.starttime<obs_time and st[0].stats.endtime>obs_time:
                    for ch in SM_chan:
                        fname=SM_dir + file.split('D')[0]+ch+'.sac'
                        st=read(fname)
                        st=st.trim(obs_time-120,obs_time+120)
                        st=st.detrend(type='demean')
                        st=st.filter('bandpass',freqmin=fmin,freqmax=fmax)
                        arr.append(st[0].trim(obs_time-120,obs_time+120))
                        #st[0].plot(starttime=obs_time,endtime=obs_time+30)
            if len(arr)==0:
                continue
            print (eve_mag,Cat['origin_time'][ieve],dist)
            print (arr[0].stats)
            tstring='%d%02d%02d%02d0000'%(Cat['origin_time'][ieve].year,Cat['origin_time'][ieve].month,Cat['origin_time'][ieve].day,Cat['origin_time'][ieve].hour)
            for ch in BK_chan:
                fname=('%s/BK.BKS.%s.%s.sac'%(BK_dir,ch,tstring))
                st=read(fname)
                st=st.trim(obs_time-60,obs_time+60)
                st=st.detrend(type='demean')
                st=st.filter('bandpass',freqmin=fmin,freqmax=fmax)
                arr.append(st[0].trim(obs_time-60,obs_time+30))
                #st[0].plot(starttime=obs_time,endtime=obs_time+30)
            S=Stream(arr)
            S.plot(equal_scale=False)
                #if st[0].stats.channel=='BHZ' and st[0].stats.starttime<obs_time and st[0].stats.endtime>obs_time:


        ieve+=1


def read_sta_file_ncedc():

    fp=open('../aux/stations_ncedc.text')
    l=fp.readlines()
    content = [x.strip() for x in l]
    fp.close()
    sta_list=[]
    line=0
    while line < len(content):

        L=content[line]
        S={}
        
        if L[0]=='#':
            line+=1
            continue
        
        S['net']=L.split('|')[0]
        S['name']=L.split('|')[1]
        
        S['loc']=L.split('|')[2]
        S['latitude']=float(L.split('|')[4])
        S['longitude']=float(L.split('|')[5])
        S['elevation']=float(L.split('|')[6])
        S['channel']=[L.split('|')[3]]
        S['unit']=[L.split('|')[13]]
        
        if "M/S**2" in S['unit']:
            S['unit']=['ACC M/S**2']
        elif "M/S" in S['unit']:
            S['unit']=['VEL M/S']
        else:
            line+=1
            continue

        S['depth']=[float(L.split('|')[7])]

        T=datetime.strptime(L.split('|')[15],'%Y-%m-%dT%H:%M:%S')
        S['starttime']=[T]
        T=datetime.strptime(L.split('|')[16],'%Y-%m-%dT%H:%M:%S')
        S['endtime']=[T]

        S['azm']=[float(L.split('|')[8])]
        S['dip']=[float(L.split('|')[9])]

        j=line+1
        while j<len(content) and content[j].split('|')[1]==S['name']:
            S['channel'].append(content[j].split('|')[3])
            T=datetime.strptime(content[j].split('|')[15],'%Y-%m-%dT%H:%M:%S')
            S['starttime'].append(T)
            T=datetime.strptime(content[j].split('|')[16],'%Y-%m-%dT%H:%M:%S')
            S['endtime'].append(T)
            S['azm'].append(float(content[j].split('|')[8]))
            S['dip'].append(float(content[j].split('|')[9]))
            u=content[j].split('|')[13]
            if "S**2" in u:
                S['unit'].append('ACC M/S**2')
            elif "S" in u:
                S['unit'].append('VEL M/S')
            else:
                j+=1
                continue
            j+=1
        
        sta_list.append(S)
        try:
            print (line,j,S['name'],content[j].split('|')[1])
        except:
            break
        line=j
    
    fp=open('../aux/sta_list_ncedc_2000-2017.pyc','wb')
    pickle.dump(sta_list,fp)
    fp.close()


def rotate(st,t0,sec_length):
    
    
    from obspy import Trace,Stream
    from obspy.geodetics.base import gps2dist_azimuth
   
    nstart=int((t0-st[0].stats.starttime)*st[0].stats.sampling_rate)
    nend=nstart+int(sec_length*st[0].stats.sampling_rate)


    if nstart<0:
        print (nstart,nend)
        return -1

    for s in st:
        if 'x' in s.stats.channel or 'X' in s.stats.channel:
            trX=np.copy(s.data[nstart:nend])
        elif 'y' in s.stats.channel or 'Y' in s.stats.channel:
            trY=np.copy(s.data[nstart:nend])

        
    B=np.reshape(np.hstack((np.repeat(trX,1),np.repeat(trY,1))),(2,len(trX)*1))

    #print (B.shape,B[:,0],B[:,1],trX[0],trX[1])
    #return
    tmp_max=0
    for deg in range(0,180,1):
        rad=np.deg2rad(deg)
        A=np.asarray(((np.cos(rad),-np.sin(rad)),(np.sin(rad),np.cos(rad))))
        C=np.dot(A,B)
        
        if np.max(np.abs(C[0,:]))>tmp_max:
            tmp_max=np.max(np.abs(C[0,:]))
            inc_rot=deg
            R=np.copy(C[0,:])
            T=np.copy(C[1,:])
        #print (D[0][0],D[0][1],D[1][0],D[1][1])

    new_stream=Stream()
    if tmp_max==0:
        return -1
    print (R.shape,T.shape,inc_rot)
    for comp in ['R','T']:
        trace=Trace()
        trace.stats=st[0].stats.copy()
        trace.stats.starttime=t0
        #trace.stats=st[0].stats
        
        if 'R' in comp :
            trace.data=R
            trace.stats.channel='xxR'
        else:
            trace.data=T
            trace.stats.channel='xxT'

        new_stream.append(trace)
    
    trace=Trace()
    trace.data=np.copy(st[2].data)
    trace.stats=st[0].stats.copy()
    trace.trim(t0,t0+sec_length-st[0].stats.delta)
    trace.stats.channel='xxZ'

    new_stream.append(trace)
    #new_stream.trim(t0,t0+sec_length)
    print (new_stream[2].stats.starttime,new_stream[2].stats.endtime)

    lat_sta=new_stream[0].stats.location.latitude
    lon_sta=new_stream[0].stats.location.longitude
    lat_eve=new_stream[0].stats.evelat
    lon_eve=new_stream[0].stats.evelon
    tmp=gps2dist_azimuth(lat_eve,lon_eve,lat_sta,lon_sta)
    
    for j in range(0,3):
        print (new_stream[j].stats.starttime,new_stream[j].stats.endtime)
        new_stream[j].stats['back_azimuth']=tmp[2]
        new_stream[j].stats['longtitude']=lon_sta
        new_stream[j].stats['latitude']=lat_sta
        new_stream[j].stats.pop('location',None)
    
    #new_stream.rotate('RT->NE')
    return new_stream

def get_trig_blast():

    from obspy import Trace,Stream
    from obspy.core.utcdatetime import UTCDateTime
    
    hd5="/home/ainbal/Data/myshakeMeta_before_20161205.h5"
    dir_data="/data/sensordata/output"
    fp=open('../aux/query_blast.txt')
    Cat={}
    Cat['origin_time']=[[]] 
    Cat['lat']=[[]]
    Cat['lon']=[[]]
    Cat['depth']=[[]]
    Cat['mag']=[[]]
    Cat['eve_id']=[[]]

    for line in fp:
        tt=datetime.strptime(line.split()[0]+line.split()[1].split('.')[0],'%Y/%m/%d%H:%M:%S')+\
                timedelta(microseconds=int(line.split()[1].split('.')[1])*1e4)
        
        Cat['origin_time'][0]=tt
        Cat['lat'][0]=(float(line.split()[2]))
        Cat['lon'][0]=(float(line.split()[3]))
        Cat['depth'][0]=(float(line.split()[4]))
        Cat['mag'][0]=(float(line.split()[5]))
        Cat['eve_id'][0]=(int(line.split()[12]))

        tt=tt-timedelta(seconds=60)
        qbeg='%4d%02d%02d %02d:%02d:%02d'%(tt.year,tt.month,tt.day,tt.hour,tt.minute,tt.second)
        tmp=tt+timedelta(seconds=300)
        qend='%4d%02d%02d %02d:%02d:%02d'%(tmp.year,tmp.month,tmp.day,tmp.hour,tmp.minute,tmp.second)

        query = 'index>\"%s \" & index<\"%s \"'%(qbeg,qend)
    
        hdf_trig = read_hdf(hd5, key = 'Trigger', where = query)
        hdf_trig = hdf_trig[['deviceid','tt', 'latitude', 'longitude', 'tf','ti']]
        
        lat_org=float(line.split()[2])
        lon_org=float(line.split()[3])
        lon0=great_circle(distance=-1e4, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
        lon1=great_circle(distance=1e4, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
        lat0=great_circle(distance=-1e4, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
        lat1=great_circle(distance=1e4, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
        
        hdf_trig = hdf_trig[(hdf_trig['latitude']>lat0) & (hdf_trig['latitude']<lat1) & (hdf_trig['longitude']>lon0) &\
                                                (hdf_trig['longitude']<lon1)]
        
        if len(hdf_trig)<1:
            continue

        Obj={}
        Obj['eve_files']={}
        Obj['eve_files']['files']=[[]]
        Obj['eve_files']['eve']=[[]]
        Obj['eve_files']['eve'][0]=({key : Cat[key][0] for key in Cat.keys()})


        for trig_time,dev,tt,ti in zip(hdf_trig.index,hdf_trig['deviceid'],hdf_trig['tt'],hdf_trig['ti']):
            path=build_path(dir_data,dev,tt)
            if path==-1:
                continue
            Obj['eve_files']['files'][0].append(path)
       
        out=read_file_from_eve(Obj,0)
        out=read_data_from_Traces(out)
        if len(out[0])==0:
            continue
        fname='tmp/%d_M%.1f.pyc'%(Cat['eve_id'][0],Cat['mag'][0])
        t0=UTCDateTime(Cat['origin_time'][0])-10

        O=[]
        for cx,cy,cz in zip(out[0],out[1],out[2]):
            o=Stream([cx,cy,cz])
            r=rotate(o,t0,40)
            if r==-1:
                continue
            O.append(r)
        
        if len(O):
            fout=open(fname,'wb')
            pickle.dump(O,fout)

    fp.close()

def get_stations_bay_area():
    
    import matplotlib as mpl
    import matplotlib.cm as cm
    import pickle
    
    lat=[]
    lon=[]
    inc=86400
    hd5="/home/ainbal/Data/myshakeMeta_before_20161205.h5"
    lon0=-123.898751872
    lon1=-119.338974772
    lat0=35.4559243504
    lat1=39.1369625317
    Id=[]

    for month in range(3,4):


        if (month==12):
            ndays=31
        else:
            ndays=(datetime(2016,month+1,1)-datetime(2016,month,1)).days
        
        nhr=ndays*24
        nsec=nhr*3600
        
        for sec in range(0,nsec,inc):
            tt0=datetime(2016,month,1)+timedelta(seconds=sec)
            tt1=tt0+timedelta(seconds=inc)
            
            print (tt0,tt1)       
            qbeg='%4d%02d%02d %02d:%02d:%02d'%(tt0.year,tt0.month,tt0.day,tt0.hour,tt0.minute,tt0.second)
            qend='%4d%02d%02d %02d:%02d:%02d'%(tt1.year,tt1.month,tt1.day,tt1.hour,tt1.minute,tt1.second)
            query = 'index>\"%s \" & index<\"%s \"'%(qbeg,qend)
            hdf_beat = read_hdf(hd5, key = 'Heartbeat', where = query)


            hdf_beat = hdf_beat[(hdf_beat['latitude']>lat0) & (hdf_beat['latitude']<lat1) & (hdf_beat['longitude']>lon0) &
                        (hdf_beat['longitude']<lon1)]
            hdf_beat = hdf_beat[['deviceId','createdOn', 'latitude', 'longitude']]
            hdf_beat = hdf_beat.drop_duplicates(['deviceId'])
            Id=Id+list(set([i for i in hdf_beat['deviceId']]))
            print (len(Id))
        
        tmax=tt1

    Id=list(set(Id))
    #print len(Id)
    fp=open('Id.pyc','wb')
    pickle.dump(Id,fp)
    fp.close()
    inc=3600
    count=0
    for month in range(3,4):
        if (month==12):
            ndays=31
        else:
            ndays=(datetime(2016,month+1,1)-datetime(2016,month,1)).days
        
        nhr=ndays*24
        nsec=nhr*3600

        for sec in range(0,nsec,inc):
            tt0=datetime(2016,month,1)+timedelta(seconds=sec)
            tt1=tt0+timedelta(seconds=inc)
            
            if tt1>tmax:
                break

            qbeg='%4d%02d%02d %02d:%02d:%02d'%(tt0.year,tt0.month,tt0.day,tt0.hour,tt0.minute,tt0.second)
            qend='%4d%02d%02d %02d:%02d:%02d'%(tt1.year,tt1.month,tt1.day,tt1.hour,tt1.minute,tt1.second)
            query = 'index>\"%s \" & index<\"%s \"'%(qbeg,qend)
            hdf_beat = read_hdf(hd5, key = 'Heartbeat', where = query)


            hdf_beat = hdf_beat[(hdf_beat['latitude']>lat0) & (hdf_beat['latitude']<lat1) & (hdf_beat['longitude']>lon0) &
                        (hdf_beat['longitude']<lon1)]
            hdf_beat = hdf_beat[['deviceId','createdOn', 'latitude', 'longitude']]
            hdf_beat = hdf_beat.drop_duplicates(['deviceId'])

           
            out=[]
            print (tt0,tt1)
            for i in range(0,len(hdf_beat)):
                indx=Id.index(hdf_beat['deviceId'][i])
                out.append([hdf_beat['latitude'][i],hdf_beat['longitude'][i],indx])
            
            out=np.asarray(out)
            day=tt0.day
            hour=tt0.hour
            if os.path.isdir('tmp/%02d/%02d'%(month,day))==0:
                os.mkdir('tmp/%02d/%02d'%(month,day))
            outfile='tmp/%02d/%02d/%02d.bin'%(month,day,hour)
            fp=open(outfile,'wb')
            out.tofile(outfile)
            fp.close()
    
    
            

def get_trig_eew(win):

    def con_line(line):
        Str=line.split('|')[5].split('Z')[0].split('.')
        if len(Str)==1:
            return Str[0]+'.000'
        else:
            return Str[0]+'.'+'%d'%(int(Str[1])*1000)


    trig_file=open('../aux/EEWHistory.txt','r')
    txt=trig_file.read().split('\n')
    trig_file.close()
    fp=open('trig_hb.txt','w')

    hd5="/home/ainbal/Data/myshakeMeta_before_20161205.h5"

    iline=0

    while iline<len(txt):
        if iline<10 or txt[iline][0]=='+':
            iline+=1
            continue

        ##NOTE
        
        if int(txt[iline].split('|')[1])<2459 or int(txt[iline].split('|')[1])>2463:
            iline+=1
            continue

        
        try:
            Str=con_line(txt[iline])
        except:
            iline+=1
            continue
        tt=datetime.strptime(Str,' %Y-%m-%dT%H:%M:%S.%f')
        eve_lat=float(txt[iline].split('|')[3])
        eve_lon=float(txt[iline].split('|')[4])
        eve_mag=float(txt[iline].split('|')[10])
        dist_c=dist_crit(eve_mag,1)

        if tt<datetime(2016,2,15):
            iline+=1
            continue
        
        jline=iline+1
        try:
            Str=con_line(txt[jline])
        except:
            iline+=1
            continue

        tt_tmp=datetime.strptime(Str,' %Y-%m-%dT%H:%M:%S.%f')
        dt_sec=((tt_tmp-tt).days*86400)+((tt_tmp-tt).seconds)
        while dt_sec<300 and jline<len(txt):
            jline+=1
            try:
                Str=con_line(txt[jline])
            except:
                jline+=1
                continue
            tt_tmp=datetime.strptime(Str,' %Y-%m-%dT%H:%M:%S.%f')
            dt_sec=(tt_tmp-tt).days*86400+(tt_tmp-tt).seconds

        if jline==len(txt):
            break
        tt=tt-timedelta(seconds=60)
        qbeg='%4d%02d%02d %02d:%02d:%02d'%(tt.year,tt.month,tt.day,tt.hour,tt.minute,tt.second)
        tmp=tt+timedelta(seconds=120)
        qend='%4d%02d%02d %02d:%02d:%02d'%(tmp.year,tmp.month,tmp.day,tmp.hour,tmp.minute,tmp.second)
        

        query = 'index>\"%s \" & index<\"%s \"'%(qbeg,qend)
        
        #print (query,iline,jline,dt_sec,tt,tt_tmp)
        hdf_trig = read_hdf(hd5, key = 'Trigger', where = query)
        hdf_trig = hdf_trig[['deviceid','tt', 'latitude', 'longitude', 'tf']]
        hdf_trig = hdf_trig[(hdf_trig['tf']==4)]
        
        trig_count=0

        for lon,lat in zip(hdf_trig['longitude'],hdf_trig['latitude']):
            dist=vincenty((eve_lat,eve_lon),(lat,lon)).kilometers
            if dist<=dist_c:
                trig_count+=1
                #print (lat,lon)

        #print (query,trig_count)
        
        trig_id = [dev for dev in  hdf_trig['deviceid']]
        
        tmp=tt-timedelta(seconds=3600)
        qend=qbeg
        qbeg='%4d%02d%02d %02d:%02d:%02d'%(tmp.year,tmp.month,tmp.day,tmp.hour,tmp.minute,tmp.second)
        #query_small = 'index>\"%s \" & index<\"%s \"'%(qbeg,qend)
        query_small='index>"20160223 02:18:51 " & index<"20160223 03:18:51 "'

        #query_small='index> 2016-02-23T02:18:53 & index< 2016-02-23T03:18:53'
        hdf_beat = read_hdf(hd5, key = 'Heartbeat', where = query_small)
        print (hdf_beat.index[0],hdf_beat.index[len(hdf_beat.index)-1],len(hdf_beat))
        #hdf_beat = hdf_beat.drop_duplicates(['deviceId'])
        hb_count=0

        for dev,lon,lat in zip(hdf_beat['deviceId'],hdf_beat['longitude'],hdf_beat['latitude']):
            dist=vincenty((eve_lat,eve_lon),(lat,lon)).kilometers
            if dist<=dist_c:
                hb_count+=1
        
        if hb_count==0:
            iline=jline
            continue
        fp.write('%s %f %f %.1f %d %d %lf\n'%(tt,eve_lat,eve_lon,eve_mag,trig_count,hb_count,trig_count/hb_count*100.0))
        #print (tt,trig_count,hb_count,trig_count/hb_count*100.0,iline)
        print (query_small,hb_count)
        iline+=1



        
        #iline=jline
        

    fp.close()    
    return
    
    qbeg=20160101
    qend=20170131
    lat_min=win['lat'][0]
    lat_max=win['lat'][1]
    lon_min=win['lon'][0]
    lon_max=win['lon'][1]

    
    j=0
    while j<len(df_trig.index):
        time=df_trig.index[j]
        lat_org=df_trig['latitude'][j]
        lon_org=df_trig['longitude'][j]
        
        
        for k in range(j+1,len(df_trig.index)):
            
            if ((df_trig.index[k]-time).days*86400+(df_trig.index[k]-time).seconds) > 30:
                break

        if k==j:
            break

        trig_id=[df_trig['deviceid'][j1] for j1 in range(j,k)]
        
        print (j,time,len(trig_id))

        t0=time-timedelta(seconds=7200)
        qbeg='%d%02d%02d %02d:%02d:%02d'%(t0.year,t0.month,t0.day,t0.hour,t0.minute,t0.second)
        qend='%d%02d%02d %02d:%02d:%02d'%(time.year,time.month,time.day,time.hour,time.minute,t0.second)
        query_small = 'index>\"%s\" & index<\"%s\"'%(qbeg,qend)
        
        lon0=great_circle(distance=-1e4, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
        lon1=great_circle(distance=1e4, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
        lat0=great_circle(distance=-1e4, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
        lat1=great_circle(distance=1e4, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
        
        hdf_beat = read_hdf(hd5, key = 'Heartbeat', where = query_small)
        hdf_beat = hdf_beat[(hdf_beat['latitude']>lat0) & (hdf_beat['latitude']<lat1) & (hdf_beat['longitude']>lon0) &
                           (hdf_beat['longitude']<lon1)]
        
        hdf_beat = hdf_beat[['deviceId', 'latitude', 'longitude']]
        hdf_beat = hdf_beat.drop_duplicates(['deviceId'])
        count=0

        for h in hdf_beat['deviceId']:
            if h in trig_id:
                print (h)
                continue
            count+=1
        
        if count>0:
            print (k-j,count)
        #if count==0:
        #    print (qbeg,qend,lon0,lon1,lat0,lat1)
        #    print (hdf_beat['deviceId'],trig_id)
        #    break
        j=k





def get_noise_tmp_coverage(win):
    
    from scipy.signal import butter, filtfilt
    from scipy.integrate import cumtrapz
    from obspy.core.utcdatetime import UTCDateTime

    def butter_bandpass(lowcut, highcut, fs, order=4):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype='band')
        return b, a


    #def butter_bandpass_filter(data, b, highcut, fs, order=5):
        #b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        #y = lfilter(b, a, data)
        #return y
    
    
    def Make_stream(data,tstart,ttrig,fmin,fmax,iproc,nsec_noise):
        
        #data are constructed:
            
   #                         -> f0
   #                 -> acc. -> .
   #                         -> .
   #                         -> fn
   #                         -> f0
   # comp.[x,y,z]    -> vel. -> .
   #                         -> .
   #                         -> fn
   #                         -> f0
   #                 -> dis. -> .
   #                         -> .
   #                         -> fn
            
            
        #NOTE
        #from obspy import Trace,Stream
        #from obspy.core.trace import Stats
        
        tmp=data.copy()
        ttrig=ttrig.to_datetime64()
        tstart=tstart.to_datetime64()
        
        if tstart+np.timedelta64(nsec_noise,'s')>=ttrig-np.timedelta64(1,'s'):
            return(0,0)
        
        nend=nsec_noise*fs-1
        g2ms=9.80655
        
        if len(data)/4/25.0 < 10.0:
            return (0,0)

        try:
            data=np.reshape(data,(4,int(len(data)/4)))[:,0:nend]
            S=np.sum(data[0:3],axis=1)
            if len(np.nonzero(S)[0])==0 or data.shape[-1]/25.0<10.0:
                return (0,0)
        except:
            return (0,0)
        
        #TIME ARRAY
        T=np.asarray([tstart+np.timedelta64(int(t*1000),'ms') for t in data[3]],dtype=np.datetime64)
        #CHECK SAMPLES RELATIVE TO TRIGGER
        nttrig=np.where(T-ttrig<=np.timedelta64(-10,'s'))[0]
        if len(nttrig)==0:
            return (0,0)

        #CHECK SAMPLES NOT MISSING
        indx=np.where(np.diff(T[nttrig]) > np.timedelta64(50,'ms'))[0]
        
        if len(indx)==0:
            indx_good=nttrig
       
        else:
            indx_good=range(0,indx[0]+1)
        
        #MIN. DURATION IS 10 S.
        if len(indx_good)/fs<10.0:
            return(0,0)
        
        data=data[:,indx_good]
        data[0:3,:]*=g2ms
        data[0:3,:]-=np.repeat(np.mean(data[0:3],axis=1),data.shape[1]).reshape(3,data.shape[1])
        
        x=np.tile(np.copy(data[3]),3).reshape((3,data.shape[1]))
        
        ifreq=0
        #DUPLICATE X,Y,Z BEFORE INTEGRATION
        data=np.tile(data[0:3],len(UNIT)).reshape(3,len(UNIT),data.shape[1])
        #INTEGRATE X,Y,Z ACC. TO GET VEL.
        data[:,1,1:]=cumtrapz(data[:,1,:],x,axis=1)
        #ASSIGN FIRST ELEMENT THE VALUE OF THE SECOND ELEMENT
        data[:,1,0]=data[:,1,1]
        #INTEGRATE X,Y,Z VEL. TO GET DISP.
        data[:,2,2:]=cumtrapz(data[:,1,1:],x[:,1:],axis=1)
        #ASSIGN FIRST TWO ELEMENTS THE VALUE OF THE THIRD ELEMENT
        data[:,2,0:2]=np.repeat(data[:,2,2],2).reshape(3,2)
        mean=np.mean(data,axis=-1)
        data-=np.repeat(mean,data.shape[-1]).reshape(data.shape)
        #print (np.mean(data,axis=-1))
        
        if SPEC_FLAG==1:
            
            O=np.zeros((len(COMP),len(UNIT),len(fmin)))

            amp=np.abs(np.fft.rfft(data))
            Freq=np.fft.rfftfreq(data.shape[-1], d=0.04)
            
            #plt.subplot(3,1,1)
            #plt.loglog(Freq,amp[0,0,:])
            #plt.ylim(1e-2,10)
            #plt.subplot(3,1,2)
            #plt.loglog(Freq,amp[0,1,:])
            #plt.ylim(1e-2,10)
            #plt.subplot(3,1,3)
            #plt.plot(x[0],data[0,1,:])
            #plt.show()
            
            rows=np.array([0,1,2],dtype=np.intp)

            for f0,f1 in zip(fmin,fmax):
                indx=np.where((Freq>f0) & (Freq<f1))
                

                if len(indx[0])==0:
                    continue
                
                col=indx[0]        
                
                chX=np.mean(amp[0,rows[:,np.newaxis],col],axis=-1)
                chY=np.mean(amp[1,rows[:,np.newaxis],col],axis=-1)
                chZ=np.mean(amp[2,rows[:,np.newaxis],col],axis=-1)
                
                if chX[0]>chY[0]:
                    O[0,:,ifreq]=chX
                else:
                    O[0,:,ifreq]=chY
                
                O[1,:,ifreq]=chZ
                ifreq+=1
            
            return (O,data.shape[-1])
         

        elif CC_FLAG==1:
            ifreq=0 
            for f0,f1 in zip(fmin,fmax):
                if f0!=3.2:
                    ifreq+=1
                    continue

                filtered=filtfilt(cB[ifreq],cA[ifreq],data)
                CC_out=np.zeros(3)
                for k in range(3):
                    CC=np.correlate(filtered[0,k]/np.std(filtered[0,k]),filtered[1,k]/np.std(filtered[1,k]),mode='valid')/len(filtered[0,k])
                    CCmax=np.max(CC)
                    CCmin=np.min(CC)
                    if np.abs(CCmax)>np.abs(CCmin):
                        CC_out[k]=CCmax
                    else:
                        CC_out[k]=CCmin
                
                return (CC_out,data.shape[-1])
        
        O=np.empty((len(COMP),len(UNIT),len(fmin),data.shape[-1]))

        for b,a in zip(cB,cA):
            
            #filter
            filtered=filtfilt(b,a,data)
            
            
            #duplicate x,y,z before integration
            #filtered=np.tile(filtered,len(UNIT)).reshape(3,len(UNIT),data.shape[1])
            
            #integrate x,y,z acc. to get vel.
            #filtered[:,1,1:]=cumtrapz(filtered[:,1,:],x,axis=1)
            
            #assign first element the value of the second element
            #filtered[:,1,0]=filtered[:,1,1]
            
            #integrate x,y,z vel. to get disp.
            #filtered[:,2,2:]=cumtrapz(filtered[:,1,1:],x[:,1:],axis=1)
            
            #assign first two elements the value of the third element
            #filtered[:,2,0:2]=np.repeat(filtered[:,2,2],2).reshape(3,2)
            #filtered=np.power(filtered,2)
            
            #HORIZONTALS
            O[0,:,ifreq]=np.sum(np.power(filtered[0:2],2),axis=0)
            #VERTICALS
            O[1,:,ifreq]=np.power(filtered[2],2)

            ifreq+=1
        
        
        
        return (O,data.shape[-1])

        print (ttrig,tstart)
        
        data=tmp.copy()

        arr_T=data[3]
        st=Stream()

        
        data=np.reshape(data,(4,int(len(data)/4)))
        ttrig=UTCDateTime(ttrig)
        tstart=UTCDateTime(tstart)
        arr_T=data[3]
        st=Stream()

        for icomp,comp in enumerate(['x','y','z']):
            
            if len(data[icomp])==0:
                return (0,0)

            tr=Trace()
            stats=Stats()
            stats.network='BM'
            stats.sampling_rate=25
            stats.channel=comp
            stats.ttrig=ttrig
            stats.starttime=tstart
            stats.npts=len(data[icomp])
            #NOTE
            tr.data=data[icomp]*g2ms
            tr.stats=stats
            st.append(tr)
        
        st.trim(tstart,tstart+nsec_noise-1.0/25.0)
        st=st.detrend(type='demean')
        
        N=np.sum([t.stats.npts for t in st])
        data_length=st[0].stats.npts

        if N==0:
            return (0,0)
        
        for t in st:
            if t.stats.npts-np.count_nonzero(t.data)==t.stats.npts:
                return (0,0)
        
        O={}

        for unit in UNIT:
            O[unit]={}
            for comp in COMP:
                O[unit][comp]=[]
                for f0,f1 in zip(fmin,fmax):
                    Str=st.copy().filter('bandpass',freqmin=f0,freqmax=f1)
                    
                    if unit=='vel':
                        Str=Str.integrate()
                    
                    if unit=='disp':
                        Str=Str.integrate().integrate()
                    
                    if comp=='hor':
                        O[unit][comp].append([Str[0].data**2+Str[1].data**2])
                    else:
                        O[unit][comp].append([Str[2].data**2])

        
        return (O,data_length)


    def check_on_land(clat,clon):

        dist_min=1e9
        for iseg,seg in enumerate(coast):
            for n,p in enumerate(seg):
                dist=(vincenty((clat,clon),(p[0],p[1])).kilometers)
                if dist<dist_min:
                    dist_min=dist
                    seg_min=iseg
                    pmin=p
        
        if clon<pmin[1]:
            return 1
        
        return 0
     

    def trig_area_par(nwin_p_proc,offset,q,iproc,hd5,win,fmin,fmax):

        import struct
       
        offset_iproc=iproc*nbin_day*nbin_hr*ncomp*nfreq
        
        out=[]
        for j in range(nbin_day):
            out.append([])
            for k in range(nbin_hr):
                out[j].append([])


        hdf_trig_local=hdf_trig
        nsec_noise=50
        count=0
        ntrace_max=10
        tot=len(hdf_trig_local)
        prec_10=int(0.1*tot)
        fp=open('test_data.bin','rb')

        t0=time.time()
        flag=0
        
        
        for i in range(offset,offset+nwin_p_proc):
            
            nb_offset=hdf_trig_local['nbytes_cum'].iloc[i]
            nb_read=hdf_trig_local['nbytes_stream'].iloc[i]

            t=(hdf_trig_local['tstart'].iloc[i].weekday(),hdf_trig_local['tstart'].iloc[i].hour)
            
            iday=t[0]
            ihour=t[1]
            fp.seek(0)
            fp.seek(nb_offset)
            data=fp.read(nb_read)
            data=np.asarray(struct.unpack(int(len(data)/4)*'f',data[:]))
           
            #seperate night- and day-time
            if nightDay==1: 
                Tlocal=hdf_trig_local['tstart'].iloc[i].tz_localize('UTC').tz_convert('US/Pacific')
                iday=0
                #nighttime
                if Tlocal.hour>=20 or Tlocal.hour<6:
                    ihour=0
                #daytime
                else:
                    ihour=1

            #stream0,data_length=Make_stream(data,hdf_trig_local['tstart'].iloc[i],hdf_trig_local['ttrig'].iloc[i],fmin,fmax,iproc,nsec_noise)
            stream0,data_length=Make_stream(data,hdf_trig_local['tstart'].iloc[i],hdf_trig_local['ttrig'].iloc[i],fmin,fmax,iproc,nsec_noise)
            
            if data_length==0:
                continue
           
            arr_dayHour[iproc*(nbin_day*nbin_hr)+iday*nbin_hr+ihour]+=data_length
            latlon_cell[iproc*(ave+1)*3+(i-offset)*3]=hdf_trig_local['latitude'].iloc[i]
            latlon_cell[iproc*(ave+1)*3+(i-offset)*3+1]=hdf_trig_local['longitude'].iloc[i]
            latlon_cell[iproc*(ave+1)*3+(i-offset)*3+2]=hdf_trig_local['tstart'].iloc[i].tz_localize('UTC').tz_convert('US/Pacific').hour
            out[iday][ihour].append(stream0)
            flag=1
            

        for iday in range(nbin_day):
            for ihour in range(nbin_hr):
                
                if len(out[iday][ihour])==0:
                    continue
                
                for icomp in range(ncomp):
                        for ifreq in range(nfreq):
                        
                            indx=(iday*nbin_hr + ihour ) * ncomp*nfreq + icomp*nfreq + ifreq
                            if SPEC_FLAG==0 and CC_FLAG==0:
                                rms_dayHour_acc[offset_iproc+indx]=np.sqrt(np.mean(np.hstack([o[icomp][0][ifreq] for o in out[iday][ihour]]),axis=-1))
                                rms_dayHour_vel[offset_iproc+indx]=np.sqrt(np.mean(np.hstack([o[icomp][1][ifreq] for o in out[iday][ihour]]),axis=-1))
                                rms_dayHour_disp[offset_iproc+indx]=np.sqrt(np.mean(np.hstack([o[icomp][2][ifreq] for o in out[iday][ihour]]),axis=-1))
                            elif CC_FLAG==1 and fmin[ifreq]==3.2 and icomp==0:
                                rms_dayHour_acc[offset_iproc+indx]=np.mean(np.hstack([o[0] for o in out[iday][ihour]]),axis=-1)
                                rms_dayHour_vel[offset_iproc+indx]=np.mean(np.hstack([o[1] for o in out[iday][ihour]]),axis=-1)
                                rms_dayHour_disp[offset_iproc+indx]=np.mean(np.hstack([o[2] for o in out[iday][ihour]]),axis=-1)
                            elif SPEC_FLAG:
                                rms_dayHour_acc[offset_iproc+indx]=np.mean(np.hstack([o[icomp][0][ifreq] for o in out[iday][ihour]]))
                                rms_dayHour_vel[offset_iproc+indx]=np.mean(np.hstack([o[icomp][1][ifreq] for o in out[iday][ihour]]))
                                rms_dayHour_disp[offset_iproc+indx]=np.mean(np.hstack([o[icomp][2][ifreq] for o in out[iday][ihour]]))

                            #if ifreq==len(fmin)-1:
                            #    print\
                            #    (iproc,COMP[icomp],fmin[ifreq],rms_dayHour_acc[offset_iproc+indx],\
                            #     rms_dayHour_vel[offset_iproc+indx],rms_dayHour_disp[offset_iproc+indx],\
                            #     len(np.hstack([o[icomp][0][ifreq] for o in out[iday][ihour]])))
        
        #for j in range(nbin_day*nbin_hr*ncomp*nfreq):
        #    if len(sub_l_acc[j])>0:
        #        
        #        tmp=np.median(np.hstack(sub_l_acc[j]))
        #        #tmp=np.sqrt(np.median(tmp[tmp>0]))
        #        rms_dayHour_acc[offset_iproc+j]=tmp
                
        #        tmp=np.median(np.hstack(sub_l_vel[j]))
        #        #tmp=np.sqrt(np.median(tmp[tmp>0]))
        #        rms_dayHour_vel[offset_iproc+j]=tmp
                
        #        tmp=np.median(np.hstack(sub_l_disp[j]))
        #        #tmp=np.sqrt(np.median(tmp[tmp>0]))
        #        rms_dayHour_disp[offset_iproc+j]=tmp
        #    else:
        #        rms_dayHour_acc[offset_iproc+j]=0
        #        rms_dayHour_vel[offset_iproc+j]=0
        #        rms_dayHour_disp[offset_iproc+j]=0
        
        print ('finished reading ',nwin_p_proc,'traces','in ',time.time()-t0,'sec.')
        return
    
    import matplotlib.path as mplPath
    import itertools

    dir_data="/data/sensordata/output"
    coast_file='LA_coast.txt'
    fp=open(coast_file,'r')
    P=[]
    for line in fp:
        if line[0]=='#' or line[0]=='>':
            continue
        
        
        clat=float(line.split()[1])
        clon=float(line.split()[0])
        if (clat<33.6 and clon<-118):
            continue
        #print (clon,clat)
        P.append([clat,clon])
    
    P.append([33.35,-117.4])
    P.append([34.5,-117.4])      
    Qhull=mplPath.Path(np.array(P))
    
    Output={}
    search_rad=[0,5,10,20,50]
    mag=np.asarray([1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5])
    nproc=8
    UNIT=['acc','vel','disp']
    
    
    time_win=[datetime(2016,2,20),datetime.now()]
    qbeg='%d%02d%02d %02d:%02d:00'%(time_win[0].year,time_win[0].month,time_win[0].day,time_win[0].hour,time_win[0].minute)
    qend='%d%02d%02d %02d:%02d:00'%(time_win[1].year,time_win[1].month,time_win[1].day,time_win[1].hour,time_win[1].minute)
    nsec=(time_win[1]-time_win[0]).days*86400
    manager = Manager()
    lon_org=win['lon'][0]
    lat_org=win['lat'][0]
    xmin=0
    xmax=(vincenty((lat_org,lon_org),(lat_org,win['lon'][1])).kilometers)
    ymin=0
    ymax=(vincenty((lat_org,lon_org),(win['lat'][1],lon_org)).kilometers)
    DX=2
    reg_name=win['reg']
    hd5="Mdata_%s.h5"%(reg_name)
    if os.path.isdir(reg_name)==0:
        os.mkdir(reg_name)


    #hd5="/home/u2/ainbal/tmp/Data/myshakeMeta.h5"

    nbinx_km=int((xmax-xmin)/DX)+1
    nbiny_km=int((ymax-ymin)/DX)+1
    
    if nightDay==0:
        nbin_day=7
        nbin_hr=24
    else:
        nbin_day=1
        nbin_hr=2
   
    #xstart=46
    #ystart=26
    ystart=yend=0
    xstart=xend=0
    sub_reg=0
    if sub_reg==1:
        fcoord=open('coord_LA_sub.lonlat','w')
        fp_coord_cell=open('cell_coord_LA_sub','wb')
        if DX==2:
            xstart=0
            xend=nbinx_km
            ystart=10
            yend=nbiny_km
        if DX==4:
            xstart=10
            xend=30
            ystart=10
            yend=40

    else:
        fcoord=open('coord_LA.lonlat','w')
        fp_coord_cell=open('cell_coord_LA','wb')

    fmin=[0.1,0.2,0.4,0.8,1.6,3.2,6.4]
    fmax=[0.2,0.4,0.8,1.6,3.2,6.4,12]
    COMP=['hor','ver']


    for c in COMP:
        DIR="%s/noise_comp_%s"%(reg_name,c)
        if os.path.isdir(DIR)==0:
            os.mkdir(DIR)

    ncomp=len(COMP)
    fs=25
    cB=[]
    cA=[]
    Be=[]

    for f0,f1 in zip(fmin,fmax):
        b,a=butter_bandpass(f0,f1,fs)
        cB.append(b)
        cA.append(a)
        Be.append((f0+f1)/2)

    Be=np.tile(np.asarray(Be),len(UNIT)).reshape(len(UNIT),len(fmin))
    SPEC_FLAG=0
    CC_FLAG=1

    nfreq=len(fmin)
    
    for iy in range(0,nbiny_km):
        
        if sub_reg and (iy<ystart or iy>yend):
            continue
    
        for ix in range(0,nbinx_km):
            
            if sub_reg and (ix<xstart or ix>xend):
                continue
            
            lon0=great_circle(distance=ix*DX*1e3, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
            lon1=great_circle(distance=(ix+1)*DX*1e3 , azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
            lat0=great_circle(distance=iy*DX*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
            lat1=great_circle(distance=(iy+1)*DX*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
            clon=great_circle(distance=DX*(ix+0.5)*1e3, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
            clat=great_circle(distance=DX*(iy+0.5)*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
            
            if (clat<33.6 and clon<-118) or Qhull.contains_point((clat,clon))==0:
                continue
            
            fcoord.write("%f %f\n"%(clon,clat))
            
            win_search={'lon':[lon0,lon1], 'lat':[lat0,lat1]}
            query_small = 'tstart>=\"%s\" & tstart<\"%s\" & latitude>=\"%s\" & latitude<\"%s\" & longitude>=\"%s\" & longitude<\"%s\"'%\
                    (qbeg,qend,win_search['lat'][0],win_search['lat'][1],win_search['lon'][0],win_search['lon'][1])
            
            hdf_trig = read_hdf(hd5, 'df',where = query_small)
            

            ave=int(len(hdf_trig)/nproc)
            extra=len(hdf_trig)%nproc
            offset=0
            
            arr_dayHour=Array(c_double , np.zeros(nproc*nbin_day*nbin_hr) , lock=multiprocessing.Lock())
            rms_dayHour_acc=Array(c_double , np.zeros(nproc*nbin_day*nbin_hr*ncomp*nfreq) , lock=multiprocessing.Lock())
            rms_dayHour_vel=Array(c_double , np.zeros(nproc*nbin_day*nbin_hr*ncomp*nfreq) , lock=multiprocessing.Lock())
            rms_dayHour_disp=Array(c_double , np.zeros(nproc*nbin_day*nbin_hr*ncomp*nfreq) , lock=multiprocessing.Lock())
            latlon_cell=Array(c_double , np.ones(nproc*(ave+1)*3)*-1 , lock=multiprocessing.Lock())
            
            #print (iy,ix,len(hdf_trig))
            if len(hdf_trig)>0:
                
                jobs=[]
                for iproc in range(0,nproc):
                    if iproc<extra:
                        nwin_p_proc=ave+1
                    else:
                        nwin_p_proc=ave
                     
                    q=multiprocessing.Queue()
                    p=multiprocessing.Process(target=trig_area_par,args=(nwin_p_proc,offset,q,iproc,hd5,win_search,fmin,fmax))
                    jobs.append([p,q])
                    p.start()
                    #offset+=timedelta(seconds=nwin_p_proc)
                    offset+=nwin_p_proc
        
                for ijob,j in enumerate(jobs):
                    j[0].join()
                    #if ijob>0:
                    #    arr_dayHour[0:nbin_day*nbin_hr]=np.add(arr_dayHour[0:nbin_day*nbin_hr],arr_dayHour[ijob*nbin_day*nbin_hr:(ijob+1)*nbin_day*nbin_hr])
            
            arr_dayHour=np.asarray(arr_dayHour)
            rms_dayHour_acc=np.asarray(rms_dayHour_acc)
            rms_dayHour_vel=np.asarray(rms_dayHour_vel)
            rms_dayHour_disp=np.asarray(rms_dayHour_disp)

            latlon_cell=np.array(latlon_cell,dtype=np.float)
            latlon_cell=latlon_cell[np.where(latlon_cell!=-1)]
            
            if len(latlon_cell)>0:
                latlon_cell=np.reshape(latlon_cell,(int(len(latlon_cell)/3),3))
                fp_coord_cell.write(b">%lf %lf, %d %d\n"%(clat,clon,iy,ix))
                np.savetxt(fp_coord_cell,latlon_cell,fmt='%lf')
            
            print ('y:',iy,"/",nbiny_km,'x:',ix,"/",nbinx_km,len(np.nonzero(arr_dayHour)[0]))
            
            for ic,co in enumerate(COMP):
                if CC_FLAG and ic!=0:
                    break

                for day in range(nbin_day):
                    for hr in range(nbin_hr):
                        for ifreq,freq in enumerate(zip(fmin,fmax)):

                            if CC_FLAG==1 and freq[0]!=3.2:
                                continue

                            if freq[0]!=freq[1]:
                                if SPEC_FLAG:
                                    Dir='%s/noise_comp_%s/DX%dkm_f%.2f-%.2f_SPEC'%(reg_name,co,DX,freq[0],freq[1])
                                elif CC_FLAG:
                                    Dir='%s/noise_comp_%s/DX%dkm_f%.2f-%.2f_CC'%(reg_name,co,DX,freq[0],freq[1])
                                else:
                                    Dir='%s/noise_comp_%s/DX%dkm_f%.2f-%.2f'%(reg_name,co,DX,freq[0],freq[1])
                            else:
                                Dir='%s/noise_comp_%s/DX%dkm_nofilt'%(reg_name,co,DX)
                            
                            if os.path.isdir(Dir)==0:
                                os.mkdir(Dir)
                            
                            tmp=[]
                            rms_acc=[]
                            rms_vel=[]
                            rms_disp=[]
                            w=[]
                            
                            if CC_FLAG:
                                if len(rms_dayHour_acc[rms_dayHour_acc>0])==0:
                                    break

                                ave_acc=np.mean(rms_dayHour_acc[rms_dayHour_acc>0])
                                ave_vel=np.mean(rms_dayHour_vel[rms_dayHour_acc>0])
                                ave_disp=np.mean(rms_dayHour_disp[rms_dayHour_acc>0])
                                fp=open('%s/%02d_%02d_f%.2f-%.2f.xyz'%(Dir,day,hr,freq[0],freq[1]),'a')
                                fp.write('%f %f %.8g %.8g %.8g\n'%(clon,clat,ave_acc,ave_vel,ave_disp))
                                continue


                            for iproc in range(nproc):
                                rms_acc.append(rms_dayHour_acc[iproc*nbin_day*nbin_hr*ncomp*nfreq + day*nbin_hr*ncomp*nfreq + hr*ncomp*nfreq + ic*nfreq + ifreq])
                                rms_vel.append(rms_dayHour_vel[iproc*nbin_day*nbin_hr*ncomp*nfreq + day*nbin_hr*ncomp*nfreq + hr*ncomp*nfreq + ic*nfreq + ifreq])
                                rms_disp.append(rms_dayHour_disp[iproc*nbin_day*nbin_hr*ncomp*nfreq + day*nbin_hr*ncomp*nfreq + hr*ncomp*nfreq + ic*nfreq + ifreq])
                                w.append(arr_dayHour[iproc*nbin_day*nbin_hr+day*nbin_hr+hr])

                            tmp=list(itertools.chain.from_iterable(tmp))
                            tmp=np.ravel(np.asarray(tmp))

                            w=np.asarray(w)
                            if np.sum(w)==0:
                                ave_acc=ave_vel=ave_disp=0
                            
                            else:
                                ave_acc=np.average(rms_acc,weights=w)
                                ave_vel=np.average(rms_vel,weights=w)
                                ave_disp=np.average(rms_disp,weights=w)
                                if co=='hor' and hr==1:
                                    print (co,hr,ifreq,ave_acc,ave_vel,ave_disp)
                                
                            if freq[0]!=freq[1]:
                                if nightDay==0:
                                    fp=open('%s/%02d_%02d_f%.2f-%.2f.xyz'%(Dir,day,hr,freq[0],freq[1]),'a')
                                else:
                                    fp=open('%s/%02d_%02d_f%.2f-%.2f.xyz'%(Dir,day,hr,freq[0],freq[1]),'a')
                            else:
                                fp=open('%s/%02d_%02d.xyz'%(Dir,day,hr),'a')
                            
                            if ave_acc>0:
                                fp.write('%f %f %.8g %.8g %.8g %d\n'%(clon,clat,ave_acc,ave_vel,ave_disp,np.sum(w)))
                            else:
                                fp.write('%f %f 0.0 0.0 0.0 0\n'%(clon,clat))

                            fp.close()

                     
            del arr_dayHour,rms_dayHour_acc,rms_dayHour_vel,rms_dayHour_disp,latlon_cell
        
    fcoord.close()
    
    return     

    #return (np.reshape(np.asarray(arr_dayHour[0:nbin_day*nbin_hr]),(nbin_day,nbin_hr)),L)

def extract_noise(win):
    
    lon_org=win['lon'][0]
    lat_org=win['lat'][0]
    xmin=0
    xmax=(vincenty((lat_org,lon_org),(lat_org,win['lon'][1])).kilometers)
    ymin=0
    ymax=(vincenty((lat_org,lon_org),(win['lat'][1],lon_org)).kilometers)
    DX=2

    nbinx_km=int((xmax-xmin)/DX)+1
    nbiny_km=int((ymax-ymin)/DX)+1
    nbin_day=7
    nbin_hr=24
   
    print (nbinx_km,nbiny_km)
    Mat=np.zeros((nbiny_km,nbinx_km,nbin_day,nbin_hr))
    coord=[]
    
    sign=0
    for iy in range(0,nbiny_km):
        for ix in range(0,nbinx_km):

            if os.path.isfile('../noise/y_%d_x_%d.pyc'%(iy,ix)) is False:
                clon=great_circle(distance=DX*(ix+0.5)*1e3, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
                clat=great_circle(distance=DX*(iy+0.5)*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
                coord.append([clon,clat])
                #for day in range(nbin_day):
                #    for hour in range(nbin_hr):
                #        Mat[iy][ix][day][hour]=np.nan
                continue

            Obj=pickle.load(open('../noise/y_%d_x_%d.pyc'%(iy,ix),'rb'))
            coord.append([Obj['clon'],Obj['clat']])
            
            for data in Obj['acc']:
                day=data['day']
                hr=data['hour']
                comp=data['comp']
                
                if comp!='x' or data['data'].shape[0]==0 or len(data['data'].shape)==1:
                    #Mat[iy][ix][day][hr]=np.nan
                    continue
                #mean=np.mean(data['data'])
                #print (day,hr,comp,data['data'].shape,len(data['data'].flatten()))
                try:
                    Mat[iy][ix][day][hr]=np.sqrt(np.mean(np.power(data['data'].flatten(),2)))
                except:
                    print ('prob. reading data')
            
    
    for day in range(nbin_day):
        for hour in range(nbin_hr):
            
            fp=open('../noise/day_hr/%02d_%02d.xyz'%(day,hour),'w')
            #fp.write('#nx=%d ny=%d range=%f %f %f %f\n'%(nbinx_km,nbiny_km,np.min(lon),np.max(lon),np.min(lat),np.max(lat)))
            for iy in range(0,nbiny_km):
                for ix in range(0,nbinx_km):
                    fp.write('%f %f %f\n'%(coord[iy*nbinx_km+ix][0],coord[iy*nbinx_km+ix][1],Mat[iy][ix][day][hour]))
            fp.close()
    
    
def MakeInv(obj):

    ttrig=pd.to_datetime(obj['triggerTimer'],unit='ms')
    tstart=pd.to_datetime(obj['data']['ts'][0],unit='ms')
    tend=pd.to_datetime(obj['data']['ts'][len(obj['data']['ts'])-1],unit='ms')

    sta = obspy.core.inventory.Station(
        latitude=obj['location']['latitude'],
        longitude=obj['location']['longitude'],
        elevation=obj['location']['altitude'],
        code=obj['deviceId'].replace('-',''),
        start_date=tstart,
        end_date=tend,
        creation_date=obj['createdOn'],
        site=obspy.station.util.Site(''),
        channels=[])
    

    for comp in ['x','y','z']:
        if comp in obj['data'].keys():
            sta.channels.append(
                Channel(code=comp,
                        location_code="",
                        latitude=obj['location']['latitude'],
                        longitude=obj['location']['longitude'],
                        elevation=obj['location']['altitude'],
                        depth=0,
                        comments=[obspy.core.inventory.util.Comment(id=0,value='accuracy=%s'%(obj['location']['accuracy']))],
                        start_date=tstart,
                        end_date=tend))

    #net = Network(code="Berkeley-MyShake",stations=[],description="smartphone array")
    #inv = Inventory(networks=[],source="MyShake01")
    #net.stations.append(sta)
    #inv.networks.append(net)
    #inv.write("station.xml", format="stationxml", validate=True)
    return sta

def MakeData(obj):
   
    npts_use=60*25-1
    
    npts=len(obj['data']['ts'])


    ttrig=pd.to_datetime(obj['triggerTimer'],unit='ms')
    tstart=pd.to_datetime(obj['data']['ts'][0],unit='ms')
    tend=pd.to_datetime(obj['data']['ts'][-1],unit='ms')
    longitude=obj['location']['longitude']
    latitude=obj['location']['latitude']
    elevation=obj['location']['altitude']
    accuracy=obj['location']['accuracy']
    station=obj['deviceId']
    createdOn=obj['createdOn']
    
    npts=len(obj['data']['ts'])
    
    if npts<=npts_use:
        npts_use=npts-1
    
    if pd.to_datetime(obj['data']['ts'][npts_use])>ttrig:
        return []

        
    Datax=obj['data']['x']
    Datay=obj['data']['y']
    Dataz=obj['data']['z']
    arr_T=(obj['data']['ts']-obj['data']['ts'][0])*1e-3

    Mdata={'ttrig':ttrig,'tstart':tstart,'tend':tend,'longitude':longitude,
               'latitude':latitude,'elevation':elevation,'accuracy':accuracy,
           'station':station,'createdOn':createdOn}
    
    Data=[Datax[0:npts_use],Datay[0:npts_use],Dataz[0:npts_use],arr_T[0:npts_use]]

    return (Mdata,Data)

    

def MakeStream(obj):
    
        
    ttrig=obspy.UTCDateTime(pd.to_datetime(obj['triggerTimer'],unit='ms'))
    tstart=obspy.UTCDateTime(pd.to_datetime(obj['data']['ts'][0],unit='ms'))
    tend=obspy.UTCDateTime(pd.to_datetime(obj['data']['ts'][len(obj['data']['ts'])-1],unit='ms'))
    #NOTE
    #arr_T=np.asarray([obspy.UTCDateTime(pd.to_datetime(t,unit='ms'))-tstart for t in obj['data']['ts']])
    arr_T=(obj['data']['ts']-obj['data']['ts'][0])*1e-3
    st=obspy.Stream()
    for comp in ['x','y','z','ts']:
        tr=obspy.Trace()
        stats=Stats()
        stats.network='BM'
        stats.sampling_rate=25
        stats.longitude=obj['location']['longitude']
        stats.latitude=obj['location']['latitude']
        stats.elevation=obj['location']['altitude']
        stats.accuracy=obj['location']['accuracy']
        stats.station=obj['deviceId'].replace('-','')
        stats.channel=comp
        stats.ttrig=ttrig
        stats.createdOn=obj['createdOn']
        stats.starttime=tstart
        stats.npts=len(obj['data']['ts'])
        #stats.ts=arr_T
        if comp=='ts':
            tr.data=arr_T
        else:
            tr.data=obj['data'][comp]
        tr.stats=stats
        st.append(tr)

    return st


def check_in_hd5(obj):
  
    tstart=obspy.UTCDateTime(pd.to_datetime(obj['data']['ts'][0],unit='ms'))
    devId=obj['deviceId'].replace('-','')
    try: 
        with pyasdf.ASDFDataSet("new_file.h5", compression="gzip-3", mode='r') as HD5:
            for s in HD5.ifilter(HD5.q.station==devId,HD5.q.starttime==tstart):
                return 1

            return 0
    except:
        return 0


def read_write_hd5(win):
  
    import struct
    import mmap

    def trig_area_par2hd5(t0,t1,q,iproc,hd5):
        
        #time_win=[t0,t0+timedelta(seconds=nwin_p_proc)]
        time_win=[t0,t1]
        
        qbeg='%d%02d%02d %02d:%02d:00'%(time_win[0].year,time_win[0].month,time_win[0].day,time_win[0].hour,time_win[0].minute)
        qend='%d%02d%02d %02d:%02d:00'%(time_win[1].year,time_win[1].month,time_win[1].day,time_win[1].hour,time_win[1].minute)
        DX=20
        print (iproc,time_win,qbeg,qend)
        
        #query_small = 'index>\"%s\" & index<\"%s\"'%(qbeg,qend)
        
        for iy in range(0,nbiny_km):
            for ix in range(0,nbinx_km):

                lon0=great_circle(distance=ix*DX*1e3, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
                lon1=great_circle(distance=(ix+1)*DX*1e3 , azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
                lat0=great_circle(distance=iy*DX*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
                lat1=great_circle(distance=(iy+1)*DX*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
                query_small = 'index>=\"%s\" & index<\"%s\" & latitude>=\"%s\" & latitude<\"%s\" & longitude>=\"%s\" &longitude<\"%s\"'%(qbeg,qend,lat0,lat1,lon0,lon1)
                hdf_trig = read_hdf(hd5, key = 'Trigger', where = query_small)
                #NOTE
                tot=len(hdf_trig.index)
                sub_l=[]
                count=0
                datafile='tmp_%s/ip%d/%02d_%02d.bin'%(reg,iproc,iy,ix)

                if os.path.isfile(datafile) or len(hdf_trig.tt.values)==0:
                    continue

                for item in zip_longest(hdf_trig.index, hdf_trig.tt, hdf_trig.deviceid):
                    tt=item[0]
                    tstamp=item[1]
                    dev=item[2]

                    t=(tt.weekday(),tt.hour)
                    iday=t[0]
                    ihour=t[1]
                    
                    year=tt.year
                    month=tt.month
                    day=tt.day
                    hour=tt.hour
                    path='%s/%d/%02d/%02d/%02d:00:00/%s_%s.pkl.gz'%(dir_data,year,month,day,hour,dev,tstamp)
                     
                    try:
                        inF = gzip.open(path, "rb")
                    except:
                        continue
                    
                    obj=pickle.load(inF,encoding='latin1')
                    inF.close()
                    
                    stream0=MakeData(obj)
                    if len(stream0)==0:
                        continue

                    if int(((count*100)/tot))>0 and count%100==0:
                        print (iproc,count,tot,int(((count*100)/tot)),'%',iy,ix)

                    sub_l.append([stream0[0],stream0[1]]) 
                    count+=1

                if len(sub_l)==0:
                    print ('file ',datafile,'does not exist','len sub_l=',len(sub_l))
                    continue
                
                fdata=open(datafile,'wb')
                Mdata=[]
                for obj in sub_l:
                    Data=[]
                    for k in range(0,4):
                        Data.append(obj[1][k])
                    
                    Data=np.hstack(Data)
                    p=struct.pack(len(Data)*'f',*(Data.astype('f4')))
                    fdata.write(p)
                    nbytes=len(p)
                    obj[0]['nbytes_offset']=nbytes
                    Mdata.append(obj[0])
                fdata.close()
                
                fdata=open('tmp_%s/ip%d/%02d_%02d.pyc'%(reg,iproc,iy,ix),'wb')
                pickle.dump(Mdata,fdata)
                fdata.close()
                

            
    
    hd5="/home/u2/ainbal/tmp/Data/myshakeMeta.h5"
    dir_data="/data/sensordata/output"
    poly_file='../aux/LA_metro.lonlat'
    fp=open(poly_file,'r')

    
    nproc=8
    ncomp=3
    
    time_win=[datetime(2016,2,20),datetime.now()]
    qbeg='%d%02d%02d %02d:%02d:00'%(time_win[0].year,time_win[0].month,time_win[0].day,time_win[0].hour,time_win[0].minute)
    qend='%d%02d%02d %02d:%02d:00'%(time_win[1].year,time_win[1].month,time_win[1].day,time_win[1].hour,time_win[1].minute)
    
    query_small = 'index>=\"%s\" & index<\"%s\" & latitude>=\"%s\" & latitude<\"%s\" & longitude>=\"%s\"\
            &longitude<\"%s\"'%(qbeg,qend,win['lat'][0],win['lat'][1],win['lon'][0],win['lon'][1])

    hdf_trig = read_hdf(hd5, key = 'Trigger', columns=['latitude'],where = query_small)
    time_proc=np.array_split(hdf_trig.index,nproc)
    del hdf_trig

    nsec=(time_win[1]-time_win[0]).days*86400
    ave=int(nsec/nproc)
    extra=nsec%nproc
    jobs=[]
    offset=time_win[0]
    lon_org=win['lon'][0]
    lat_org=win['lat'][0]
    xmin=0
    xmax=(vincenty((lat_org,lon_org),(lat_org,win['lon'][1])).kilometers)
    ymin=0
    ymax=(vincenty((lat_org,lon_org),(win['lat'][1],lon_org)).kilometers)
    DX=20
    reg_name=win['reg']

    nbinx_km=int((xmax-xmin)/DX)+1
    nbiny_km=int((ymax-ymin)/DX)+1

    for iproc in range(0,nproc):
        q=multiprocessing.Queue()
        p=multiprocessing.Process(target=trig_area_par2hd5,args=(time_proc[iproc][0],time_proc[iproc][-1],q,iproc,hd5))
        jobs.append([p,q])
        p.start()
        #print(offset,offset+timedelta(seconds=nwin_p_proc))
        #offset+=timedelta(seconds=nwin_p_proc)
    
    #NOTE 
    for ijob,j in enumerate(jobs):
        j[0].join()
    
    D=0
    fdata=open('test_data_%s.bin'%(reg),'wb')
    
    print ('start writing binary datafile')
    batch=1

    for iy in range(0,nbiny_km):
        for ix in range(0,nbinx_km):
            print (iy,"/",nbiny_km,ix,"/",nbinx_km)
            for iproc in range(nproc):

                try:
                    M=pickle.load(open('tmp_%s/ip%d/%02d_%02d.pyc'%(reg,iproc,iy,ix),'rb'))
                except:
                    continue
                
                if D==0:
                    D={}
                    keys=M[0].keys()
                    for key in keys:
                        D[key]=[]

                for key in keys:
                    for obj in M:
                        D[key].append(obj[key])
                
                print ('read D, start read tmp_%s/ip%d/%02d_%02d.bin'%(reg,iproc,iy,ix))
                
                fp=open('tmp_%s/ip%d/%02d_%02d.bin'%(reg,iproc,iy,ix),'rb')
                sum_c=0 
                for obj in M:
                    nbytes=obj['nbytes_offset']
                    count=int(nbytes/4)
                    
                    if nbytes==0:
                        continue

                    Data=np.fromfile(fp,dtype='f4',count=count)
                    sum_c+=nbytes
                    if len(Data)!=count:
                        print (len(Data),count)
                        exit()
                    
                    n=int(len(Data)/4)
                    try:
                        Min=np.min(Data[n*3:])
                    except:
                        print(len(Data),count,nbytes,sum_c)
                        exit()

                    if Min!=0:
                        print (np.min(Data[n*3:]))
                        exit()

                    p=struct.pack(len(Data)*'f',*(Data.astype('f4')))
                    fdata.write(p)
                del Data

    fdata.close()
    print ('closed bin. datafile')
    
    tmp=np.concatenate(([0],D['nbytes_offset']))
    nbytes_cum=np.cumsum(tmp)
    D['nbytes_cum']=nbytes_cum[0:len(nbytes_cum)-1]
    D['nbytes_stream']=D.pop('nbytes_offset')
    df = pd.DataFrame(D)
    store=pd.HDFStore('new_%s.hd5'%(reg),mode='w')
    store.append('Mdata_%s'%(reg),df,data_columns=True,index=False)
    print ('added Mdata to store')
    store.close()



if __name__ == '__main__':
    reg='LA'
    win=get_spat_win(reg)
    #read_write_hd5(win)
    #exit()
    #win=get_spat_win('LA')
    #print (win)
    nightDay=1
    get_noise_tmp_coverage(win)

    #read_write_hd5(win)

