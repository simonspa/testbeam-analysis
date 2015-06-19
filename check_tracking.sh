#!/bin/bash

logfiles=/tmp/logs

# Delete old log files:
rm -r $logfiles

# Inflate all log archives:
for file in `ls output/logs/ | grep analysis`
do 
    unzip -qq output/logs/$file -d $logfiles
    grep -e "run header" -e "GBL tracks" -e "DUT links" $logfiles/${file//zip/log}
done

