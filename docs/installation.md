
## Install with pypi, conda or docker

### PyPI

```bash
pip install trackplot
# __Note:__ We noticed some pypi mirrors are not syncing some packages we depend on, 
# therefore please try another pypi mirror once you encounter 
# `No local packages or working download links found for xxx`

pip install pybigwig hicmatrix # to enable bigWig, bigBed and hicMatrix support
```

### Conda

```bash
conda install -c bioconda -c conda-forge trackplot

# or install trackplot into an isolated environments
conda create -n trackplot -c bioconda -c conda-forge trackplot

# or install latest trackplot
git clone https://github.com/ygidtu/trackplot.git trackplot
cd trackplot
conda create -n trackplot -f environment.yaml
```

### Docker

```bash
docker pull ygidtu/trackplot
docker run --rm ygidtu/trackplot --help

# or build docker image from source
git clone https://github.com/ygidtu/trackplot.git trackplot
cd trackplot
docker build -t ygidtu/docker .
docker run --rm ygidtu/trackplot --help
```

## Install from source

download source code using git
```bash
git clone https://github.com/ygidtu/trackplot.git trackplot
cd trackplot
```

### Run as command line tools

```bash
# optional, enable bigWig, bigBed and hicMatrix support
pip install pybigwig hicmatrix

python setup.py install

trackplot --help
# or
python main.py --help
```

### Run as script

The following 3 installation methods will try to install pyBigWig and hicmatrix by default, 
once facing installation issues please check their official document.

#### 1. using python

```bash
pip install -r requirements.txt
python main.py
```

#### 2. using pipenv

```bash
pipenv install  # create virtualenv and install required packages
# optional, with `--pypi-mirror https://pypi.tuna.tsinghua.edu.cn/simple` to specify your faverate PyPi mirror
# optional, with `--skip-lock` once encounter locking issues

pipenv shell && python main.py --help    # switch to virtualenv

# or just run with pipenv
pipenv run python main.py --help
```

#### 3. using poetry

```bash
# once facing installation issues, please try to change PyPi mirror in tool.poetry.source section of pyproject.toml 
poetry install   # create virtualenv and install required packages
poetry shell  && python main.py --help   # switch to virtualenv

# or just run with poetry
poetry run python main.py --help
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
docker pull ygidtu/trackplot
docker run --rm ygidtu/trackplot --help
```

### Build docker image from source

```bash
git clone https://github.com/ygidtu/trackplot.git trackplot
cd trackplot
docker build -t ygidtu/docker .
docker run --rm ygidtu/trackplot --help
```

### Command line Usage

```bash
docker run --rm -v $PWD:$PWD --user $(id -u):$(id -g) ygidtu/trackplot --help
```

- `-v`, `--volumn`: mount the working directory to docker container
- `--user`: prevent docker read and write file using root privileges
- the rest usage please check [Command line usage](./command.md)

**Note: ** detailed command line usage please check [Command line Usage](./command.md)


## Build Web interface from source

1. **nodejs is required**
2. **Users could change the server ip and port by modify the settings.ini**

```bash
git clone https://github.com/ygidtu/trackplot.git trackplot
cd trackplot/web

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
docker pull ygidtu/trackplotweb

# -v map the current working directory into docker containers
# -p map the outer port to inner port of docker container
docker run --name trackplotweb \
  --rm -v $PWD:$PWD \
  -p 5000:5000 \
  ygidtu/trackplotweb 
```