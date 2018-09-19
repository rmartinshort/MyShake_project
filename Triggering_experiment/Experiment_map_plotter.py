#!/usr/bin/env python
#RMS 2018

#Plot the results of a triggering experiment. For each triggering test, a map, trigger counts and four histograms of delay times 
#are plotted. This is designed to give a sense of the number of phones likely to trigger when the command is sent.

import pandas as pd
import glob
import os
import numpy as np
import datetime
from dateutil import parser
import pytz
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.basemap import Basemap

# Plot style setup
from matplotlib import rcParams
rcParams["axes.labelsize"]= 24              
rcParams["font.size"]= 30
rcParams["font.weight"] = 'bold'
rcParams["axes.titlesize"]= 30              
rcParams["legend.fontsize"]= 16             
rcParams["xtick.labelsize"]= 20
rcParams["ytick.labelsize"]= 20
rcParams["legend.fontsize"]= 20
rcParams["axes.linewidth"]  = 2 
rcParams["xtick.major.size"]  = 5
rcParams["ytick.major.size"]  = 5 
rcParams["xtick.major.width"]  = 2
rcParams["ytick.major.width"]  = 2
rcParams['figure.figsize'] = (10,7)
plt.style.use("ggplot")


local_tz = pytz.timezone("US/Pacific")
utc_tz = pytz.timezone("utc")

def main():

    cwd = os.getcwd()

    #Boundary box of the region. 
    coords = [-118.7,-117.2,33.0,34.5] #For LA 

    #NOTE that we expect MongoDB files to be present in this folder
    experiment_folder = './Experiment_0912_LA_mongodb'
    os.chdir(experiment_folder)
    #Get a list of all the files containing the phones we chose to trigger
    phones_to_trigger = glob.glob("Phones_to_trigger*")

    trig_dfs = []
    trig_fnames = []
    for totrigger in phones_to_trigger:
        
        to_trigger_df = pd.read_csv(totrigger,index_col=0)
        trigger_time = totrigger.split('_')[-1].split('.')[0]
        did_trigger = glob.glob("Devices_triggered_%s*" %trigger_time)
        mongo_db_time = glob.glob("MongoDBtrigger_%s*" %trigger_time)
        try:
            fnames = [totrigger,did_trigger[0],mongo_db_time[0]]
            fnames_df = [pd.read_csv(element).reset_index(drop=True) for element in fnames]
            fnames.append(trigger_time)
            trig_fnames.append(fnames)
            trig_dfs.append(fnames_df)
        except:
            continue 

    figure_count = 1
    for i in range(len(trig_dfs)):
        
        trigger_time = int(trig_fnames[i][-1])
        trigger_time_corr = local_tz.localize(datetime.datetime.fromtimestamp(trigger_time),is_dst=None).astimezone(pytz.utc)
        
        df_trio = trig_dfs[i]
        to_trig = df_trio[0]
        triggered = df_trio[1]
        mongo = df_trio[2]
        triggered_ids = triggered['deviceId']
        
        #----------------------------------
        #Create a column that indicates whether or not a phone has been triggered
        
        triggered_ids = triggered['deviceId'].values
        all_ids = to_trig['deviceId'].values
        trigs = np.zeros(len(all_ids))
        
        for i in range(len(all_ids)):
            idval = all_ids[i]
            if idval in triggered_ids:
                trigs[i] = 1
                
        to_trig['triggered'] = trigs
            
        #----------------------------------
        
        #Note that tt is the trigger time
        mongo = mongo[mongo['deviceid'].isin(triggered_ids)]
        
        #Get the trigger delay times from the mongoDB
        mongo['tt_delay'] = mongo['tt'].apply(lambda x: 
                ((local_tz.localize(datetime.datetime.fromtimestamp(x/1000),is_dst=None).astimezone(pytz.utc)) - trigger_time_corr).total_seconds())
        
        #Get the 'database injected' delay times from the mongoDB
        mongo['ti_delay'] = mongo['ti'].apply(lambda x: 
                ((local_tz.localize(datetime.datetime.fromtimestamp(x/1000),is_dst=None).astimezone(pytz.utc)) - trigger_time_corr).total_seconds())
        
        #Get the trigger delay from the phones that triggered
        triggered['trigger_delay_triggerTimer'] = triggered['triggerTimer'].apply(lambda x: 
                      (utc_tz.localize(parser.parse(x),is_dst=None) - trigger_time_corr).total_seconds())
        
        #Get the 'created on' delay from the phones that triggered
        triggered['trigger_delay_createdOn'] = triggered['createdOn'].apply(lambda x: 
                      (utc_tz.localize(parser.parse(x),is_dst=None) - trigger_time_corr).total_seconds())
        
        #Get all the triggered phones
        triggered_phones = triggered[triggered['deviceId'].isin(to_trig['deviceId'])]
        triggered_phones = triggered_phones[triggered_phones['deviceId'].isin(triggered_phones['deviceId'].unique())]
        
        #Get all delays delays associated with the phones
        mongodb_delays_tt = []
        mongodb_delays_ti = []
        triggered_lats = []
        triggered_lons = []
        phone_delays_tT = []
        phone_delays_cO = []
        
        for i in range(len(triggered_phones)):
            phoneid = triggered_phones['deviceId'].values[i]
            mongo_delay_tt = mongo[mongo['deviceid']==phoneid]['tt_delay'].values
            mongo_delay_ti = mongo[mongo['deviceid']==phoneid]['ti_delay'].values
            phone_delay_tT = triggered_phones['trigger_delay_triggerTimer'].values[i]
            phone_delay_cO = triggered_phones['trigger_delay_createdOn'].values[i]
            phonelon = triggered_phones['longitude'].values[i]
            phonelat = triggered_phones['latitude'].values[i]
            
            #print(mongo_delay_tt,phone_delay_tT)
            #print(mongo_delay_ti,phone_delay_cO)
            try:
                mongodb_delays_tt.append(mongo_delay_tt[0])
                mongodb_delays_ti.append(mongo_delay_ti[0])
                triggered_lats.append(phonelat)
                triggered_lons.append(phonelon)
                phone_delays_tT.append(phone_delay_tT)
                phone_delays_cO.append(phone_delay_cO)
            except:
                print("No delay time!")
        
        generate_figure(coords, trigger_time_corr, to_trig, mongodb_delays_tt, 
            mongodb_delays_ti, phone_delays_tT, phone_delays_cO, figure_count)
        figure_count += 1


def generate_figure(coords,trigger_time_corr,phones_to_trigger,Md_tt,Md_ti,Ph_tt,Ph_co,counter):
    
    '''
    Function to generate figure showing the phones that have been triggered and various trigger delay times
    '''
    
    lonmin = coords[0]
    lonmax = coords[1]
    latmin = coords[2]
    latmax = coords[3]
    
    figure = plt.figure(figsize=(20,20))
    ax1 = figure.add_subplot(321)
    
    proportion_triggered = np.sum(phones_to_trigger['triggered'])/len(phones_to_trigger)
    total_phones = len(phones_to_trigger)
    sns.countplot(x='appVersion',data=phones_to_trigger,hue='triggered',ax=ax1)
    ax1.set_xlabel('App version')
    ax1.set_ylabel('Count')
    ax1.legend(loc='upper left',title='Triggered')
    ax1.grid()
    ax1.set_title('Proportion triggered: %.2f' %proportion_triggered)
    
    ax2 = figure.add_subplot(322)
    
    triggered_phones = phones_to_trigger[phones_to_trigger['triggered']==1]
    untriggered_phones = phones_to_trigger[phones_to_trigger['triggered']==0]

    map = Basemap(llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, 
                  lat_0=latmin+(latmax-latmin)/2, lon_0=lonmin+(lonmax-lonmin)/2, resolution = 'i',ax=ax2)
    
    map.fillcontinents(color='#cc9955', lake_color='aqua', zorder = 0)
    map.drawcoastlines(color = '0.15')
    xtrig, ytrig = map(triggered_phones['longitude'].values,triggered_phones['latitude'].values)
    xutrig, yutrig = map(untriggered_phones['longitude'].values,untriggered_phones['latitude'].values)
    map.scatter(xutrig,yutrig,8,marker='o',color='r')
    map.scatter(xtrig,ytrig,8,marker='o',color='g')
    ax2.set_title("Map of phones: Green = triggered, Red = did not trigger")
    
    ax3 = figure.add_subplot(323)

    Md_ti = np.array(Md_ti)
    Md_tt = np.array(Md_tt)
    Ph_tt = np.array(Ph_tt)

    Md_tt[Md_tt>20] = 0
    Ph_tt[Ph_tt>20] = 0
    Md_ti[Md_ti>20] = 0
    
    median_delay = np.median(Md_tt)
    ax3.hist(Md_tt,bins=20)
    ax3.axvline(median_delay,color='g',linestyle='--')
    ax3.set_xlim(0,10)
    ax3.set_ylabel('Number of phones')
    ax3.set_xlabel('MDB triggertime delay time [s]')
    ax3.grid()
    ax3.set_title('Delay times: MDB tt. Median: %g s' %median_delay)
    
    ax4 = figure.add_subplot(324)
    
    median_delay = np.median(Ph_tt)
    ax4.hist(Ph_tt,bins=20)
    ax4.axvline(median_delay,color='g',linestyle='--')
    ax4.set_xlim(0,10)
    ax4.set_ylabel('Number of phones')
    ax4.set_xlabel('SDB triggertime delay time [s]')
    ax4.grid()
    ax4.set_title('Delay times: SDB tt. Median: %g s' %median_delay)
    
    ax5 = figure.add_subplot(325)
    
    median_delay = np.median(Md_ti)
    ax5.hist(Md_ti,bins=20)
    ax5.axvline(median_delay,color='g',linestyle='--')
    ax5.set_xlim(0,10)
    ax5.set_ylabel('Number of phones')
    ax5.set_xlabel('MDB triggerinject delay time [s]')
    ax5.grid()
    ax5.set_title('Delay times: MDB ti. Median: %g s' %median_delay)
    
    ax6 = figure.add_subplot(326)
    
    median_delay = np.median(Ph_co)
    ax6.hist(Ph_co,bins=50)
    ax6.axvline(median_delay,color='g',linestyle='--')
    ax6.set_ylabel('Number of phones')
    ax6.set_xlabel('SBD CreatedOn delay time [s]')
    ax6.grid()
    ax6.set_title('Delay times: SDB createdOn. Median: %g s' %median_delay)
    
    figure.suptitle('Experiment from trigger at %s. Total phones: %i' %(trigger_time_corr,total_phones) )
    figname = "Triggering_experiment_%s.png" %counter
    plt.savefig(figname,dpi=400)

    print("Done figure generation!")

if __name__ == '__main__':

    main()
