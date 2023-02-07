# Monitoring

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
