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

Download [OpenJDK like Zulu](https://www.azul.com/downloads/?version=java-8-lts&os=ubuntu&architecture=x86-64-bit&package=jdk)
    
    sudo apt install ./zulu8.68.0.19-ca-jdk8.0.362-linux_amd64.deb
