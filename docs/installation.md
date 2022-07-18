
## Install from source

download source code using git
```bash
git clone from https://github.com/ygidtu/sashimi
cd sashimi
```

### Run as command line tools

```bash
python setup.py install
sashimi --help
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

### Note:
If there is any problem with installation of `cairocffi`

- please install the requirements according to the [Official Documentation of cairocffi](https://cairocffi.readthedocs.io/en/stable/overview.html)
- please try to use another backend, including `Agg`, `Pdf`, `Svg`, etc..., instead of `Cairo` by default.

---

## Using docker image

### Pull from docker hub
```bash
docker pull ygidtu/sashimi
docker run --rm ygidtu/sashimi --help
```

### Build docker image from source
```bash
git clone from https://github.com/ygidtu/sashimi
cd sashimi
docker build -t ygidtu/docker .
docker run --rm ygidtu/sashimi --help
```