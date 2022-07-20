
# Install from source

download source code using git
```bash
git clone from https://github.com/ygidtu/sashimi
cd sashimi
```

## Run as command line tools

```bash
python setup.py install
sashimi --help
```

## Run as script
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

## Note:
If there is any problem with installation of `cairocffi`

- please install the requirements according to the [Official Documentation of cairocffi](https://cairocffi.readthedocs.io/en/stable/overview.html)
- please try to use another backend, including `Agg`, `Pdf`, `Svg`, etc..., instead of `Cairo` by default.

---

# Using docker image

For users who wish to running this program on `Windows` or `macOS (Apple Silicon)`, we strongly recommend docker image.

## Pull from docker hub
```bash
docker pull ygidtu/sashimi
docker run --rm ygidtu/sashimi --help
```

## Build docker image from source
```bash
git clone from https://github.com/ygidtu/sashimi
cd sashimi
docker build -t ygidtu/docker .
docker run --rm ygidtu/sashimi --help
```

## Command line Usage

```bash
docker run --rm -v $PWD:$PWD --user $(id -u):$(id -g) ygidtu/sashimi --help
```

- `-v`, `--volumn`: mount the working directory to docker container
- `--user`: prevent docker read and write file using root privileges
- the rest usage please check [Command line usage](./command_line_usage.md)

**Note: ** detailed command line usage please check [Command line Usage](./command_line_usage.md)