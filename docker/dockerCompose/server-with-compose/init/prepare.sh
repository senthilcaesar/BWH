#!/usr/bin/env sh

rm /data/index.html
echo "<h1>Welcome from Docker Compose!</h1>" >> /data/index.html
echo "<img src='https://www.docker.com/wp-content/uploads/2022/03/Moby-logo.png' />" >> /data/index.html
