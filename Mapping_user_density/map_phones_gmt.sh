#!/bin/bash

#RMS 2018
#Map of phone user locations.
#Most mapping will be done in python, but in some cased GMT maps are preferable, so we need 
#something set up to produce them 

gmt gmtset BASEMAP_TYPE_FANCY
gmt gmtset FONT_ANNOT_PRIMARY 20

phone_locations=/Users/rmartinshort/Documents/Berkeley/MyShake_project/Triggering_stats/Experiment_0907_LA/Phone_locations.dat

ps=phone_users.ps

J=M8i
R=-118.7/-117.2/33.0/34.5 #LA region 

#Run pscoast to generate the coastlines
gmt pscoast -Rd$R -J$J -X2i -BWSne -B0.4f0.4 -Df -Slightblue -P -Lf-117.5/34.2/34.2/30.0+l"Distance [km]" --PS_MEDIA=Custom_14ix14i -K > $ps


#Plot all phone locations  
echo "Plotting phone locations"
gmt psxy $phone_locations -J$J -Sc0.1 -Gred -Rd$R -Wthinnest -O -V -Bg90/g180 >> $ps

#Post-processing stuff - make pdf
gmt ps2raster $ps -P -Tf

open phone_users.pdf