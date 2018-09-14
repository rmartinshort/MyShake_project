#!/usr/bin/env python
#RMS 2018

import numpy as np
import pandas as pd
import geopy.distance
import utm
import matplotlib.pyplot as plt
import glob
import multiprocessing as mp

#For reading and processing data
import gzip
import json

#for timing
import time

#for writing binary files 
import struct
import pickle

def MakeData(json_file):

    '''
    Convert json file into data vectors
    '''

    infile = gzip.open(json_file, "rb")
    obj=json.loads(infile.read().decode('utf-8'))

    ttrig=pd.to_datetime(obj['triggerTimer'],unit='ms')
    tstart=pd.to_datetime(obj['data'][0]['ts'],unit='ms')
    tend=pd.to_datetime(obj['data'][-1]['ts'],unit='ms')
    longitude=obj['location']['longitude']
    latitude=obj['location']['latitude']
    elevation=obj['location']['altitude']
    accuracy=obj['location']['accuracy']
    station=str(obj['deviceId'])
    createdOn=pd.to_datetime(obj['createdOn'])

    #time array
    arr_T=[(t['ts']-obj['data'][0]['ts'])*1e-3 for t in obj['data']]
    #components of the the seismogram
    Datax=np.asarray([t['x'] for t in obj['data']])
    Datay=np.asarray([t['y'] for t in obj['data']])
    Dataz=np.asarray([t['z'] for t in obj['data']])

    #Note that the semple rate migth vary across phones and time. Thus we 
    #need to do some basic processing to make sure all the traces have the same
    #length and sample rate


    #print(len(Datax),len(Datay),len(Dataz),len(arr_T))
    #print ((tend-tstart).total_seconds())
    #print(np.mean(np.diff(arr_T)))

    Mdata={'ttrig':ttrig,'tstart':tstart,'tend':tend,'longitude':longitude,
               'latitude':latitude,'elevation':elevation,'accuracy':accuracy,
           'station':station,'createdOn':createdOn}
    Data=[Datax,Datay,Dataz,arr_T]

    return Mdata,Data 


def divide_work(array,nworkers):
    
    '''Split an array across nworkers. Each worker gets some grid cells to process'''
    
    subarrays = np.array_split(array,nworkers)
    
    return subarrays

def process_worker(grid_cells,all_triggers,output,pos):
    
    '''One worker processing some collection of grid cells'''

    number_of_users = []
        
    for element in grid_cells:

        lower_left = element[0]
        upper_right = element[2]
        old_lon = lower_left[0]
        old_lat = lower_left[1]
        new_lon = upper_right[0]
        new_lat = upper_right[1]
        index = element[-1]
        i = index[0]
        j = index[1]
             
        users_between = all_triggers[(all_triggers['latitude'] >= old_lat) & (all_triggers['latitude'] <= new_lat) & 
                    (all_triggers['longitude'] >= old_lon) & (all_triggers['longitude'] <= new_lon)]

        #Here we'll need to do some calculations to, for example, process the waveform data
        data_locations = users_between['tsobj'].apply(lambda x: x[37:]+'.json.gz').values

        cell_metadata_name = 'metadata_cell_%03i_%03i.pkl' %(i,j)
        cell_data_name  = 'data_cell_%03i_%03i.bin' %(i,j)
        #load the data file


        if len(data_locations) > 0:
            print("Number of files to process inside grid cell %i,%i : %i" %(i,j,len(data_locations)))

            #Write metadata and data to pickle and binary files. There will be one file
            #per grid cell that contains data 

            cell_data_file = open(cell_data_name,'ab')
            cell_metadata_file = open(cell_metadata_name,'ab')

            Mdata = []

            for f in data_locations:

                traces_Mdata,traces = MakeData(f)

                Data = np.array([])

                for k in range(4):
                    Data = np.concatenate((Data,traces[k]))

                p=struct.pack(len(Data)*'f',*(Data.astype('f4')))
                nbytes = len(p)
                traces_Mdata['nbytes_offset']=nbytes
                Mdata.append(traces_Mdata)
                cell_data_file.write(p)

            cell_data_file.close()
            pickle.dump(Mdata,cell_metadata_file)
            cell_metadata_file.close()
            del Mdata,p

        
        number_of_users.append(len(users_between))
        
    output.put((pos,number_of_users))

def Generate_grid(dim_x,dim_y,y1,x1,z1,utm_zone,Wd,Ld,theta,phi):

    '''
    Generate grid with cells that measure dim_x by dim_y
    '''

    n_x = dim_x
    n_y = dim_y
    grid_points = []

    lx = len(np.arange(n_x,Ld,n_x))
    ly = len(np.arange(n_y,Wd,n_y))

    #Fill this grid with the number of phone users in each cell
    image_grid = np.zeros([ly,lx])

    fig = plt.figure(figsize=(10,10))

    i = 0
    for y_step in np.arange(n_y,Wd,n_y):
        y_step_dir = y_step*np.sin(theta)
        
        #Generate upper and lower bound for y
        y_new = y1 + y_step_dir
        y_old = y_new - n_y 
        j = 0 
        for x_step in np.arange(n_x,Ld,n_x):
            x_step_dir = x_step*np.cos(phi)
            
            #Generate upper and lower bound for y
            x_new = x1 + x_step_dir
            x_old = x_new - n_x
            
            #Area of the grid cell in km^2
            area = (abs(y_new-y_old)*abs(x_new-x_old))/1e6
            
            (old_lat,old_lon) = utm.to_latlon(x_old,y_old,z1,utm_zone)
            (new_lat,new_lon) = utm.to_latlon(x_new,y_new,z1,utm_zone)
            
            grid_points.append([(old_lon,old_lat),(new_lon,old_lat),(new_lon,new_lat),(old_lon,new_lat), (i,j)])
            j += 1
        i += 1

    return grid_points, image_grid


def main(region_code='LA'):

    t0 = time.time()

    nprocs = 8 
    process_output = mp.Queue()

    print("We have %i nodes available" %mp.cpu_count())

    if nprocs > mp.cpu_count():
        raise ValueError("nprocs > cpu count on this machine!")

    regions = {'LA':[33,34.5,-118.7,-117.2]}

    #Set up the coordinates of the region we wish to investiage
    region = regions[region_code]
    lat_min = region[0]
    lat_max = region[1]
    lon_min = region[2]
    lon_max = region[3]

    #lower left 
    x1,y1,z1,u = utm.from_latlon(lat_min,lon_min)
    #upper left 
    x2,y2,z2,u = utm.from_latlon(lat_max,lon_min)
    #upper right 
    x3,y3,z3,u = utm.from_latlon(lat_max,lon_max)
    #lower right
    x4,y4,z4,u = utm.from_latlon(lat_min,lon_max)

    #Distance along width (latitude)
    Wd = geopy.distance.vincenty((lat_min,lon_min),(lat_max,lon_min)).m
    #Distance along length (longitude)
    Ld = geopy.distance.vincenty((lat_min,lon_min),(lat_min,lon_max)).m

    #Direction along width (latitude)
    theta = np.arctan((y2-y1)/(x2-x1))
    #Direction along length (longitude)
    phi = np.arctan((y4-y1)/(x4-x1))

    #all_data = "/Users/rmartinshort/Documents/Berkeley/MyShake_project/Triggering_stats/trigger_databases/Devices_triggered_all_since_2016.csv"
    all_data_trigger6 = "Devices_triggered_large_request.csv"

    triggers_all = all_data_trigger6
    print("Loading metadata")
    all_triggers = pd.read_csv(triggers_all,index_col=0)

    #### Determine grid points 

    grid_points,image_grid = Generate_grid(2000, 2000, y1, x1, z1, u, Wd, Ld, theta, phi)

    #### Do calcualtion in parallel on nprocs

    subarrays = divide_work(grid_points,nprocs)
    processes = [mp.Process(target=process_worker,args=(subarrays[i],all_triggers,process_output,i)) for i in range(nprocs)]

    for p in processes:
        print('starting process %s' %p)
        p.start()

    for p in processes:
        print("joining process %s" %p)
        p.join()

    results = [process_output.get() for p in processes]

    results.sort()

    #Gather the results into a single vector

    t1 = np.array(results[0][1])
    for element in results[1:]:
        e = element[1]
        t1 = np.concatenate((t1,np.array(e)))

    #Fill the image grid 

    for i in range(len(t1)):
        indi,indj = grid_points[i][-1]
        nusers = t1[i]
        if nusers == 0:
            nusers = np.nan
        image_grid[indi,indj] = nusers

    plt.figure(figsize=(10,10))
    plt.imshow(image_grid,extent=(lon_min,lon_max,lat_max,lat_min),cmap='jet',aspect='equal')
    plt.gca().invert_yaxis()
    plt.colorbar()
    plt.savefig("test.pdf",dpi=400)

    t1 = time.time()

    print('Total time: %i' %(t1-t0))


if __name__ == '__main__':

    main()
