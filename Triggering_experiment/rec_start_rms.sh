#!/bin/bash

#This script if for sending a list of phone IDs to be triggered. Another script is then used to go and 
#fetch the data that has been produced by this process

#Send a list of phones (ids in test_dev.txt) to be triggered
#test_param.txt contains the parameters to pass to each phone. The only one that we'll really want to
#change is the length of recording time


echo `pwd`
echo "sending dev. parameters"
echo " "
curl -u "myshakeResearch:lemy\$eeYoreD" -XPOST -k --data @test_param.txt https://bsl-myshake01:8162/api/message?destination=topic://myshake.earthquakeParamConfig
echo "sending dev. list"
echo " "
curl -u "myshakeResearch:lemy\$eeYoreD" -XPOST -k --data @Phones_to_trigger.txt https://bsl-myshake01:8162/api/message?destination=topic://myshake.deviceConfig
echo " " 
