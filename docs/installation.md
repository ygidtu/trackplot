
## Install from source

download source code using git
```bash
git clone from https://github.com/ygidtu/sashimi.py sashimi
cd sashimi
```

### Run as command line tools

```bash
pip install sashimipy

# or install from source
python setup.py install

sashimipy --help
```

### Run as script
1. using pipenv
```bash
pipenv install
pipenv run python main.py --help
```

2. using python
```bash
pip install -r requirements.txt
python main.py
```

** Note: **
If there is any problem with installation of `cairocffi`

- Please install the requirements according to the [Official Documentation of cairocffi](https://cairocffi.readthedocs.io/en/stable/overview.html)
- Or try to use another backend, including `Agg`, `Pdf`, `Svg`, etc..., instead of `Cairo` by default.

But:
- the `Agg`, `PDF`, etc. backends may have problems with the small protein domains, so use as appropriate.
- ![](imgs/cmd/1.svg)

---

## Using docker image

For users who wish to running this program on `Windows` or `macOS (Apple Silicon)`, we **strongly** recommend docker image.

### Pull from docker hub
```bash
docker pull ygidtu/sashimi
docker run --rm ygidtu/sashimi --help
```

### Build docker image from source
```bash
git clone from https://github.com/ygidtu/sashimi.py sashimi
cd sashimi
docker build -t ygidtu/docker .
docker run --rm ygidtu/sashimi --help
```

### Command line Usage

```bash
docker run --rm -v $PWD:$PWD --user $(id -u):$(id -g) ygidtu/sashimi --help
```

- `-v`, `--volumn`: mount the working directory to docker container
- `--user`: prevent docker read and write file using root privileges
- the rest usage please check [Command line usage](./command_line_usage.md)

**Note: ** detailed command line usage please check [Command line Usage](./command_line_usage.md)


## Build Web interface from source

1. **nodejs is required**
2. **Users could change the server ip and port by modify the settings.ini**

```bash
git clone https://github.com/ygidtu/sashimi.py sashimi
cd sashimi/web

# build the frontend static files
npm install -g vue-cli vite && npm install
vite build

# prepare the backend server
pip install fastapi pydantic jinja2 uvicorn

python server.py --help
```


### Using docker image

#### Pull from web

```bash
docker pull ygidtu/sashimiweb

export PORT=8080  # the port your are trying to exposed
docker run --rm -v $PWD:$PWD --user $(id -u):$(id -g) -p $PORT:5000 ygidtu/sashimiweb
```