# Monitoring
When you bring multiple Shiny apps to your end users, it may be interesting to track usage of the different applications over time. It may a.o. help to understand your user base and to prioritize maintenance work for the different applications

## 1) Prometheus
Download [Prometheus monitoring system](https://prometheus.io/download/)
```
tar -xvzf prometheus-2.42.0.linux-amd64.tar.gz
cp prometheus.yml prometheus-2.41.0.linux-amd64
```
To start the ShinyProxy monitoring run
```
cd prometheus-2.41.0.linux-amd64
./prometheus --web.listen-address=0.0.0.0:7070
```
Now go to your browser and open http://localhost:7070/targets, this page should contain one target called shinyproxy with the UP state.

## 2) Export machine metrics
Download [node_exporter](https://prometheus.io/download/#:~:text=94f1fa4cd28f057c4f16dd0718acfe5bf0b5dc8185177142c6f345d8799b11b4-,node_exporter,-Exporter%20for%20machine)
```
tar -xvzf node_exporter-1.5.0.linux-amd64.tar.gz
```
To start exporting metrics run
```
cd node_exporter-1.5.0.linux-amd64
./node_exporter
```
## 3) Setup Grafana 
Now that Prometheus scrapes the metrics, we need a way to visualize these metrics.
Install Grafana
```
sudo apt-get install -y adduser libfontconfig1
wget https://dl.grafana.com/enterprise/release/grafana-enterprise_9.3.6_amd64.deb
sudo dpkg -i grafana-enterprise_9.3.6_amd64.deb
```
Next, access Grafana on http://localhost:3000 and login using admin:admin. Next add Prometheus as datasource, remember that Prometheus is running at http://localhost:7070
