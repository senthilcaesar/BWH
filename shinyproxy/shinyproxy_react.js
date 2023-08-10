npm install prometheus-query

  const prom = new PrometheusDriver({
    endpoint: "http://remnrem.net:7070/",
    baseURL: "/api/v1", // default value
  });

  const q =
    'shinyproxy_absolute_apps_running{instance="localhost:9090", job="shinyproxy"}';
  prom
    .instantQuery(q)
    .then((res) => {
      const series = res.result;
      console.log("Series", series.length);
      // series.forEach((serie) => {
      //   console.log("Serie:", serie.metric.toString());
      //   console.log("Time:", serie.value.time);
      //   console.log("Value:", serie.value.value);
      // });
    })
    .catch(console.error);
