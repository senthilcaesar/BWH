cd /root/Programme/node_exporter-1.5.0.linux-amd64
./node_exporter


cd /root/Programme/prometheus-2.41.0.linux-amd64
./prometheus --web.listen-address="0.0.0.0:7070"

# https://prometheus.io/download/ ( Additonal metrics plugins )
