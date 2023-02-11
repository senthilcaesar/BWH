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
