Enter the following command to create a new builder, which we’ll call mybuilder:

```python
docker buildx create --name mybuilder --use --bootstrap
```
You can inspect a new builder by entering
```go
docker buildx inspect mybuilder
```
