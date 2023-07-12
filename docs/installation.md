
## Install with pypi, conda or docker

### Notes
>1. For **windows**, **mac(apple silicon)** and **other arm** users, due to several requirements pysam, pybigwig and hicmatrix do not support those platforms, pleas use docker image as alternative 
>2. if `segment fault` with multiple processing, please try to use docker image, or just run with `-p 1`.
>3. if `Please install pyBigWig and hicmatrix` occurs, please check the official document of 
    [pyBigWig](https://github.com/deeptools/pyBigWig) and [hicmatrix](https://github.com/deeptools/HiCMatrix) 
    to solve their requirements.

### PyPi

```bash
# optional, enable bigWig, bigBed and hicMatrix support
pip install pybigwig hicmatrix

pip install trackplot
# __Note:__ We noticed some pypi mirrors are not syncing some packages we depend on, 
# therefore please try another pypi mirror once you encounter 
# `No local packages or working download links found for xxx`
```

### AppImage (Linux x86_64 platform only)

Due to the limitation of AppImage technic itself, we only provide AppImage for linux (x86_64 platform) users.
Once you have installation issues and not familiar with docker, 
please download the AppImage file from our releases.

All the AppImage files were tested on the official pre-built GNU/Linux distributions docker images:

- Arch `appimagecrafters/tests-env:archlinux-latest`
- Fedora `appimagecrafters/tests-env:fedora-30`
- Debian `appimagecrafters/tests-env:debian-stable`
- Ubuntu `appimagecrafters/tests-env:ubuntu-bionic`
- Centos `appimagecrafters/tests-env:centos-7`

> **Note:** the AppImage will decompress all bundled files before execution, 
> therefore it will a little bit slower than command line tools and source code

```bash
# example with version v0.2.6, please using your interested version according to your needs
export VERSION=0.2.6
chmod +x trackplot-${VERSION}-x86_64.AppImage
./trackplot-${VERSION}-x86_64.AppImage --help
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

### Install from source

download source code using git
```bash
git clone https://github.com/ygidtu/trackplot.git trackplot
cd trackplot
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

