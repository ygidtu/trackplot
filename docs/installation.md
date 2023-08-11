
## Install with pypi, conda or docker

### Notes
>1. For **windows**, **mac(apple silicon)** and **other arm** users, due to several requirements pysam, pybigwig and hicmatrix do not support those platforms, pleas use docker image as alternative 
>2. if `segment fault` with multiple processing, please try to use docker image, or just run with `-p 1`.
>3. if `Please install pyBigWig and hicmatrix` occurs, please check the official document of 
    [pyBigWig](https://github.com/deeptools/pyBigWig) and [hicmatrix](https://github.com/deeptools/HiCMatrix) 
    to solve their requirements.

### PyPi

Before running this command line, please check python (>=3.8) was installed.

```bash
# optional, enable bigWig, bigBed and hicMatrix support
pip install pybigwig hicmatrix

pip install trackplot
# __Note:__ We noticed some pypi mirrors are not syncing some packages we depend on, 
# therefore please try another pypi mirror once you encounter 
# `No local packages or working download links found for xxx` 
```
---
### AppImage (Linux/WSL x86_64 platform only)

All the AppImage files were tested on the official pre-built GNU/Linux distributions docker images:
- Arch: `appimagecrafters/tests-env:archlinux-latest`
- Fedora: `appimagecrafters/tests-env:fedora-30`
- Debian: `appimagecrafters/tests-env:debian-stable`
- Ubuntu: `appimagecrafters/tests-env:ubuntu-bionic`
- Centos: `appimagecrafters/tests-env:centos-7`

>Due to the limitation of AppImage technic itself, we only provide AppImage for Linux and Windows subsystem for linux (x86_64 platform) users.
Once you have installation issues and not familiar with docker, 
please download the AppImage file from our releases.
> 
> Once the AppImage file couldn't work properly please open an issue in this repo, 
and provide us the system platform and full error messages for us to debug.
> 
> **Notes:** 
> 1. the AppImage will decompress all bundled files before execution, 
> therefore it will a little bit slower than command line tools and source code
> 2. please use absolute path instead of relative path.

```bash
# example with version v0.3.2, please using your interested version according to your needs
export VERSION=0.3.2
chmod +x trackplot-${VERSION}-x86_64.AppImage
./trackplot-${VERSION}-x86_64.AppImage --help
```

---


### Docker

> Known issue: the logging time may have several hours mismatch with your local time, due to timezone settings inner the image.

```bash
docker pull ygidtu/trackplot
docker run --rm ygidtu/trackplot --help

# or build docker image from source
git clone https://github.com/ygidtu/trackplot trackplot
cd trackplot
docker build -t ygidtu/docker .

docker run --rm -v $PWD:$PWD --user $(id -u):$(id -g) ygidtu/trackplot --help

```
- `-v`, `--volume`: mount the working directory to docker container
- `--user`: prevent docker read and write file using root privileges
- the rest usage please check [Command line usage](./command.md)

**Note: ** detailed command line usage please check [Command line Usage](./command.md)


---

### Install from source

Prior to installing the tool from the source code, users should verify their Python version (>=3.8).

```bash
python --version
# Python 3.10.8
```

#### python3 is not available 

If your Python version does not match the requirements of Trackplot, users could follow the cmd to install and more detailed instruction for different system please refer to [here](https://realpython.com/installing-python/).  

```bash
# If the default Python version does not meet the requirements of Trackplot, 
# users should install an appropriate version. 

# 1. install from source code
#   1.1 Step 1: Download the Source Code 
wget https://www.python.org/ftp/python/3.10.12/Python-3.10.12.tgz
tar -xvzf Python-3.10.12.tgz

#   1.2 Step 2: Build Python
find ./Python-3.10.12/Python -type d | xargs chmod 0755
cd Python-3.10.12
./configure --prefix=$PWD/Python-3.10.12/Python
make
make install

# 2. install trackplot
git clone https://github.com/ygidtu/trackplot trackplot

$PWD/Python-3.10.12/python -m pip install -e trackplot/

# 3. check trackplot
$PWD/Python-3.10.12/Python-3.10.12/Python/bin/trackplot --help

```

#### python3 is available

```bash
# 1. download the trackplot
git clone https://github.com/ygidtu/trackplot trackplot
cd trackplot

# 2. install the trackplot and it's requirements to python env
pip install -e .
trackplot --help  # or python main.py --help

# Note: pybigwig hicmatrix were optional, used to enable bigWig, bigBed and hicMatrix support
# Once the installation of pybigwig and hicmatrix fails, and these two formats are not necessary, 
# you still can using trackplot in the following way
pip install -r requirements.txt
python main.py --help
```

---

### Conda

First make sure your conda is properly installed.

```bash
# Check if conda has been successfully installed.
conda --version

# if conda is not installed, refer to https://conda.io/projects/conda/en/latest/user-guide/install/download.html

```

After successful installation

```bash

# install trackplot into the default conda env 
conda install -c bioconda -c conda-forge trackplot

# or install trackplot into an isolated environments
conda create -n trackplot -c bioconda -c conda-forge trackplot

# if conda is getting stuck at solving environment', please refer to https://stackoverflow.com/a/66963979

# or install latest trackplot  
git clone https://github.com/ygidtu/trackplot.git trackplot
cd trackplot
conda env update -n trackplot --file environment.yaml

# activate the trackplot environment and execute the command line tool
conda activate trackplot
trackplot --help
```

---

### for `pipenv` or `poetry` users

> Install [pipenv](https://pipenv.pypa.io/en/latest/) or [poetry](https://python-poetry.org)  

```bash
git clone https://github.com/ygidtu/trackplot
cd trackplot

# pipenv
# create virtualenv and install required packages
pipenv install
# optional, with `--pypi-mirror https://pypi.tuna.tsinghua.edu.cn/simple` to specify your faverate PyPi mirror
# optional, with `--skip-lock` once encounter locking issues

# switch to virtualenv
pipenv shell && python main.py --help

# or just run with pipenv
pipenv run python main.py --help


# poetry
# once facing installation issues, please try to change PyPi mirror in tool.poetry.source section of pyproject.toml 
# create virtualenv and install required packages
poetry install

# switch to virtualenv
poetry shell  && python main.py --help

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
### Build Web interface from source

For this tutorial, we will assume that you have already installed Conda on your system.

1.nodejs (18.14.0 LTS above)

```shell
# if conda is getting stuck at solving environment', please refer to https://stackoverflow.com/a/66963979

conda install -c conda-forge nodejs
```

Or user could download and install nodejs from [here](https://nodejs.org/en/).

2.install from source code

Before this, please make sure that Trackplot has been properly installed in the same env.


```shell
git clone https://github.com/ygidtu/trackplot trackplot
cd trackplot/web

# build the frontend static files; If npm was not found, please install nodejs(https://nodejs.org).
npm install -g vue-cli vite && npm install
vite build

# check whether the ui is successfully compiled
ls ../ui

# then install trackplot from source code and start server as follow
python main.py --start-server --host 127.0.0.1 --port 5000 --plots ./plots
```