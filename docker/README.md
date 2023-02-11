Enter the following command to create a new builder, which we’ll call mybuilder:

```python
docker buildx create --name mybuilder --use --bootstrap
```
You can inspect a new builder by entering
```python
docker buildx inspect mybuilder
```
You can also see what runtime platforms your current builder instance supports by running
```python
docker buildx inspect --bootstrap
```
Now, you’ll jumpstart your multi-architecture build with the single docker buildx command shown below:
```python
 docker buildx build --platform=linux/arm64,linux/amd64 --push --tag remnrem/moonlight-multiarch:latest -f Dockerfile .

```
The docker `buildx build` subcommand has a number of flags which determine where the final image will be stored. By default, i.e. if none of the flags are specified, the resulting image will remain captive in docker’s internal build cache. This is unlike the regular `docker build` command which stores the resulting image in the local `docker images` list.

We can check the image with the imagetools subcommand which confirms the architecture versions that are included in the image:
```python
docker buildx imagetools inspect "remnrem/moonlight"
docker pull remnrem/moonlight
docker inspect --format "{{.Architecture}}" "remnrem/moonlight"
```
