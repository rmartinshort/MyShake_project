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
import pickle
#import _pickle as pickle
import time
from itertools import count
#from numba import jit
import multiprocessing
from multiprocessing import Manager
from multiprocessing.sharedctypes import Value, Array
from ctypes import c_double
import gzip
import os
from obspy import Trace,Stream
#from obspy.core.inventory import Inventory, Network, Station, Channel, Site
#import obspy
#import h5py
import matplotlib.path as mplPath
#import pyasdf
from itertools import *
import sys


def get_spat_win(region):

    #NOTE THAT I USED Scal HERE
    if region is 'LA':
        lat_min=33
        lat_max=34.5
        lon_min=-118.7
        lon_max=-117.2

    elif region is 'Ncal':
        lat_min=35.455
        lat_max=39.136
        lon_min=-123.898
        lon_max=-119.338

    elif region is 'SFBay':
        lat_min=36.921
        lat_max=38.494
        lon_min=-123.128
        lon_max=-121.387

    else:
        print ('error: region name not listed ... exit')
        exit()
    
    spatial_win={}
    spatial_win['lat']=[lat_min,lat_max]
    spatial_win['lon']=[lon_min,lon_max]
    spatial_win['reg']=region

    return spatial_win

def MakeData(obj):
    
    ttrig=pd.to_datetime(obj['triggerTimer'],unit='ms')
    tstart=pd.to_datetime(obj['data'][0]['ts'],unit='ms')
    tend=pd.to_datetime(obj['data'][-1]['ts'],unit='ms')
    longitude=obj['location']['longitude']
    latitude=obj['location']['latitude']
    elevation=obj['location']['altitude']
    accuracy=obj['location']['accuracy']
    station=str(obj['deviceId'])
    createdOn=pd.to_datetime(obj['createdOn'])
    #arr_T=(obj['data']['ts']-obj['data']['ts'][0])*1e-3
    arr_T=[(t['ts']-obj['data'][0]['ts'])*1e-3 for t in obj['data']]
    npts=len(arr_T)
    Datax=np.asarray([t['x'] for t in obj['data']])
    Datay=np.asarray([t['y'] for t in obj['data']])
    Dataz=np.asarray([t['z'] for t in obj['data']])

    Mdata={'ttrig':ttrig,'tstart':tstart,'tend':tend,'longitude':longitude,
               'latitude':latitude,'elevation':elevation,'accuracy':accuracy,
           'station':station,'createdOn':createdOn}
    Data=[Datax,Datay,Dataz,arr_T]

    return (Mdata,Data)


def read_write_hd5(win):
  
    import struct

    def trig_area_par_read(nwin_p_proc,t0,q,iproc,hd5,win,L):
        
        import json
        time_win=[t0,t0+timedelta(seconds=nwin_p_proc)]

        qbeg='%d%02d%02d %02d:%02d:00'%(time_win[0].year,time_win[0].month,time_win[0].day,time_win[0].hour,time_win[0].minute)
        qend='%d%02d%02d %02d:%02d:00'%(time_win[1].year,time_win[1].month,time_win[1].day,time_win[1].hour,time_win[1].minute)
        store='../%s/Mdata.h5'%(win['reg'])
        print (qbeg,qend)
        
        sub_l=[]
        offset_iproc=iproc*nbin_day*nbin_hr
        for j in range(nbin_day*nbin_hr):
            sub_l.append(manager.list(L[offset_iproc+j]))
    
        query_small = 'index>\"%s\" & index<\"%s\"'%(qbeg,qend)
        hdf_trig = read_hdf(hd5, key = 'Trigger', where = query_small)
        mask = ((hdf_trig.index>=time_win[0]) & (hdf_trig.index<time_win[1]) & (hdf_trig['latitude']>win['lat'][0]) & (hdf_trig['latitude']<win['lat'][1]) & (hdf_trig['longitude']>win['lon'][0]) &
                           (hdf_trig['longitude']<win['lon'][1]))
        #NOTE
        tot=len(hdf_trig[mask].index)
        count=0
        
        
        for item in zip(hdf_trig[mask].index, hdf_trig[mask].tt, hdf_trig[mask].deviceid):
            
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
            path='%s/%d/%02d/%02d/%02d:00:00/%s_%s.json.gz'%(dir_data,year,month,day,hour,dev,tstamp)
            try:
                inF = gzip.open(path, "rb")
            except:
                continue
            
            obj=json.loads(inF.read().decode('utf-8'))
            inF.close()

            #check=check_in_hd5(obj)
            
            #if check==1:
            #    continue

            #arr_dayHour[iproc*(nbin_day*nbin_hr)+iday*nbin_hr+ihour]+=1
             

            #if int(((count*100)/tot))>0 and count%100==0:
            #print (iproc,count,tot,int(((count*100)/tot)),'%',iy,ix)
            stream0=MakeData(obj)
            toffset=(iday*nbin_hr+ihour)
            sub_l[toffset].append([stream0[0],stream0[1]]) 
            count+=1

        #NOTE
        for j in range(nbin_day*nbin_hr):
            L[offset_iproc+j]=sub_l[j]
        
        del sub_l
        del hdf_trig,mask
        print (iproc,'end',tot)


    #WRITE NAME OF HD5
    hd5="myshakeMeta.h5"
    if len(hd5)==0:
        print ('erorr: please enter path to new hd5 file... exit')
        exit()
    
    dir_data="/data/sensordata/output"
    #poly_file='../aux/LA_metro.lonlat'
    #fp=open(poly_file,'r')
    
    #number of threads and components
    nproc=8
    ncomp=3
   
    #time window to process
    time_win=[datetime(2017,9,20),datetime(2017,9,23)]
    print ('working on '+'%d/%d/%d'%(time_win[0].year,time_win[0].month,time_win[0].day)+' < time win. < '+'%d/%d/%d'%(time_win[1].year,time_win[1].month,time_win[1].day))

    #time window for each thread
    nsec=(time_win[1]-time_win[0]).days*86400
    ave=int(nsec/nproc)
    extra=nsec%nproc
    manager = Manager()
    lon_org=win['lon'][0]
    lat_org=win['lat'][0]
    xmin=0
    xmax=(vincenty((lat_org,lon_org),(lat_org,win['lon'][1])).kilometers)
    ymin=0
    ymax=(vincenty((lat_org,lon_org),(win['lat'][1],lon_org)).kilometers)
    #sub-window in spatial domain
    DX=20
    reg_name=win['reg']

    nbinx_km=int((xmax-xmin)/DX)+1
    nbiny_km=int((ymax-ymin)/DX)+1
    nbin_day=7
    nbin_hr=24
    
    #ystart=25
    #xstart=45
    #print (nbinx_km,nbiny_km)
    fsize_old=0
    st=Stream() 
    #inv = Inventory(networks=[],source="MyShake01")
    #net = Network(code="BM",stations=[],description="smartphone array")
    sign_write_file=0

    #open and close files to remove old content
    fdata=open('test_data.bin','wb')
    fp_mdata=open('test_data.pyc','wb')
    fdata.close()
    fp_mdata.close()
    
    ystart=0
    xstart=0
    
    for iy in range(0,nbiny_km):
        for ix in range(0,nbinx_km):
            
            L = manager.list([[]]*nproc*nbin_day*nbin_hr)
            offset=time_win[0]
            lon0=great_circle(distance=ix*DX*1e3, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
            lon1=great_circle(distance=(ix+1)*DX*1e3 , azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
            lat0=great_circle(distance=iy*DX*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
            lat1=great_circle(distance=(iy+1)*DX*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
            clon=great_circle(distance=DX*(ix+0.5)*1e3, azimuth=90, latitude=lat_org, longitude=lon_org)['longitude']
            clat=great_circle(distance=DX*(iy+0.5)*1e3, azimuth=0, latitude=lat_org, longitude=lon_org)['latitude']
            
            #spatial window to process
            win_search={'lon':[lon0,lon1], 'lat':[lat0,lat1], 'reg': reg_name}
            print (iy,ix,win,lon0,lon1,lat0,lat1)
            #arr_dayHour=Array(c_double , nproc*nbin_day*nbin_hr , lock=multiprocessing.Lock())
            
            jobs=[]
            for iproc in range(0,nproc):
                if iproc<extra:
                    nwin_p_proc=ave+1
                else:
                    nwin_p_proc=ave
                

                q=multiprocessing.Queue()
                p=multiprocessing.Process(target=trig_area_par_read,args=(nwin_p_proc,offset,q,iproc,hd5,win_search,L))
                jobs.append([p,q])
                p.start()
                #print(offset,offset+timedelta(seconds=nwin_p_proc))
                offset+=timedelta(seconds=nwin_p_proc)
            
            #NOTE 
            for ijob,j in enumerate(jobs):
                j[0].join()
          
            print ('start loop on L')
            starttime=time.time()
            #NOTE
            Mdata=[]
            p=[]
            fdata=open('test_data.bin','ab')
            fp_mdata=open('test_data.pyc','ab')
            
            for icount,item in enumerate(L):
                for obj in item:
                    Data=np.asarray([])
                    nbytes=0
                    for k in range(0,4):
                        Data=np.concatenate((Data,obj[1][k]))
                    #NOTE
                    p=struct.pack(len(Data)*'f',*(Data.astype('f4')))
                    fdata.write(p)
                    nbytes=len(p)
                    obj[0]['nbytes_offset']=nbytes
                    Mdata.append(obj[0])
            
            print ('end loop on L time=%s sec.'%(time.time()-starttime))
            print (len(Mdata))
            
            fdata.close()
            if len(Mdata)>0:
                pickle.dump(Mdata,fp_mdata)
            
            fp_mdata.close()
            del L,Mdata,p


def conv_dicto2table():

    import pandas as pd
    from pandas import HDFStore, read_hdf, Timedelta
    import pickle
    import numpy as np
    import os

    objs=[]
    count=0
    fp=open('test_data.pyc','rb')
    #loop on pickled meta-data
    while 1:                           
        try:          
            objs.append(pickle.load(fp))
            count+=1
            print (count)
        except EOFError:
            break

    fp.close()

    #arrange in new dict.
    keys=objs[0][0].keys()
    D={}
    for key in keys:
        D[key]=[]
        for o in objs:
            for item in o:
                D[key].append(item[key])

                
        print (key,len(D[key]),type(D[key][0]))

    #compute cum. bytes and rename field
    tmp=np.concatenate(([0],D['nbytes_offset']))
    nbytes_cum=np.cumsum(tmp)
    D['nbytes_cum']=nbytes_cum[0:len(nbytes_cum)-1]
    D['nbytes_stream']=D.pop('nbytes_offset')

    #write hd5, enable indexing of data coloums for fast quary
    pd.set_option('io.hdf.default_format','table')
    df = pd.DataFrame(D)
    #df.to_pickle('Mdata_LA.pyc')
    store = pd.HDFStore('Mdata_%s.h5'%(reg_name),mode='w')
    store.append('df', df,data_columns=True,index=False)
    store.close()
    os.remove('test_data.pyc')


if __name__ == '__main__':
    reg_name='LA'
    win=get_spat_win(reg_name)
    read_write_hd5(win)
    conv_dicto2table()
    #read_write_hd5(win)

