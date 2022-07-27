# Web UI

## Deployment

### Build from source

#### 1. requirements

Except the sashimi, we still required `fastapi`, `uvicorn` and `pydantic`

```bash
# python requirements
pip install fastapi uvicorn pydantic

# generating web static files

cd web && yarn build

cd ..

python server.py
```

### Using docker image

#### Pull from web

```bash
docker pull ygidtu/sashimiweb

export PORT=8080  # the port your are trying to exposed
docker run --rm -v $PWD:$PWD --user $(id -u):$(id -g) -p $PORT:5000 ygidtu/sashimiweb
```

## Usage

1. Main page
![](imgs/web/main.png)
2. Configuration page
![](imgs/web/overview.png)
3. Choose reference
![](imgs/web/reference.png)
![](imgs/web/choose.png)
4. Add plot by set input file
![](imgs/web/add.png)
5. Draw the final plot
![](imgs/web/draw.png)
6. Preview image will display under the draw section
![](imgs/web/preview.png)
