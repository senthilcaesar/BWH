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
docker buildx build --push \
--platform=linux/amd64,linux/arm64 \
--tag your_docker_username/multi_arch_sample:buildx-latest .
```
