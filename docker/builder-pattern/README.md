Using the builder pattern to optimize the size of the Docker image
means that we have to maintain two `Dockerfiles` and one shell script.
We can use efficient alternative methods like multi-stage Dockerfile.

| Builder Pattern      | Multi-Stage Docker Builds |
| :---        |    :----:   |
| Need to maintain two `Dockerfiles` and a shell script | Title       |
| Need to copy the executables to the Docker host before copying them to the final Docker image | Text        |
