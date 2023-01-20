![](shinyproxy.png)

![](https://img.shields.io/badge/Platform-linux--64%20-blue.svg)
![](https://img.shields.io/badge/ShinyProxy-2.6.1%20-blue.svg)
![](https://img.shields.io/badge/Docker-20.10.22%20-blue.svg)
![](https://img.shields.io/badge/OpenJDK_Zulu-8%20-blue.svg)


# ShinyProxy
ShinyProxy is your favourite way to deploy Shiny apps in an enterprise context.
When deploying a Shiny application with ShinyProxy, the application is simply bundled 
as an R package and installed into a Docker image. Every time a user runs an application, 
a container spins up and serves the application.


# Installation

## 1. Java 8

Download [OpenJDK like Zulu](https://www.azul.com/downloads/?package=jdk)
    
    sudo apt install ./zulu8.68.0.19-ca-jdk8.0.362-linux_amd64.deb

## 2. Docker

Download and Install [Docker for ubuntu](https://docs.docker.com/engine/install/ubuntu/)

ShinyProxy needs to connect to the docker daemon to spin up the containers for the Shiny apps. 
By default ShinyProxy will do so on port 2375 of the docker host. In order to allow for connections on port 2375, 
the startup options need to be edited.
on an Ubuntu 16.04 LTS, 18.04 LTS and 20.04 LTS or a CentOS 7, RHEL 7, CentOS 8 and RHEL 8 system 
(or a similar system that uses systemd) with Docker installed from the Docker repositories, one can change the configuration using:
```
sudo systemctl edit docker
```
