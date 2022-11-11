#!/bin/bash

remote_conn='sq566@eristwo.partners.org'
remote_conn_pass='SonyXperia$9'
data_dir='/data/nsrr/datasets/cfs/polysomnography/edfs'

while read -r case
do
    sshpass -p ${remote_conn_pass} rsync -r -v --progress -e ssh ${remote_conn}:${data_dir}/${case}.edf .
done < gilla.txt
