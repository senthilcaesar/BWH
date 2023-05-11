#!/bin/bash

log=/root/Programme/logs/log_file.txt

# create log file or overrite if already present
printf "Log File - " > $log

# append date to log file
date >> $log

cd /root/Programme
rm -rf moonlight
git clone https://github.com/remnrem/moonlight.git
cd moonlight
sed -i "s/ENV MOONLIGHT_SERVER_MODE=0/ENV MOONLIGHT_SERVER_MODE=1/g" Dockerfile
docker buildx use mybuilder >> $log
docker buildx inspect mybuilder >> $log
docker buildx build --platform=linux/arm64,linux/amd64 --push --tag remnrem/moonlight:latest . &>> $log

# End time
date >> $log
