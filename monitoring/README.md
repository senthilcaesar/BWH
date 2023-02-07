# Monitoring

## 1) Prometheus
Download [Prometheus monitoring system](https://prometheus.io/download/)
```
tar -xvzf prometheus-2.42.0.linux-amd64.tar.gz
cp prometheus.yml prometheus-2.41.0.linux-amd64
```
To start the monitoring run
```
cd prometheus-2.41.0.linux-amd64
./prometheus --web.listen-address=0.0.0.0:7070
```
