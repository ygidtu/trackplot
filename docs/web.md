# Web UI

## Installation

### AppImage (Linux x86_64 users only)

Download the trackplotweb appimage file from our [release](https://github.com/ygidtu/trackplot/releases)

```bash
# example with version v0.2.6, please using your interested version according to your needs
export VERSION=0.2.6
chmod +x trackplotweb-${VERSION}-x86_64.AppImage
./trackplotweb-${VERSION}-x86_64.AppImage --help

# startup webserver
./trackplotweb-${VERSION}-x86_64.AppImage --host 127.0.0.1 --port 5000 --plots ./plots
```
    
**Note:** the `--plots` were required while using appimages

---

### Docker image (Strong recommended)

We also prepared a docker image of web server, uses could access this by following this,

```shell
docker pull ygidtu/trackplotweb

# Deploy the server
docker run --name trackplotweb \
  --rm \
  -v $PWD/data:/data \
  -v $PWD/plots/:/plots
  -p 5000:5000 \
  --data /data \
  --plots /plots \
  ygidtu/trackplotweb 
```

`-p`: public and private port for the server, default:5000(public):5000(private)
- `-v`, `--volume`: mount the working directory to docker container, for example, the `$PWD/data` could replace by the path to your directory contains all necessary data
- `--user`: prevent docker read and write file using root privileges
- the rest usage please check [Command line usage](./command.md)

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

# prepare the backend server
pip install flask

# before startup, check whether the trackplot properly installed
python -c "import trackplot; print(trackplot.__version__)"

cd ../
python server.py --help

Usage: server.py [OPTIONS]

Options:
  -h, --host TEXT     the ip address binding to
  -p, --port INTEGER  the port to listen on
  --plots PATH        the path to directory where to save the backend plot
                      data and logs, required while using appImage.
  --data PATH         the path to directory contains all necessary data files.
  --version           Show the version and exit.
  --help              Show this message and exit.
```
---

## Usage

### 1. Main page

The main page of the server, user could click the button to create the plot.

![](imgs/web/1.png)

### 2. Configuration page

At this page, user could see a progress bar at the first tip indicating the progress of the current job.

And next, user should fill in the region of interest at second tip and click the confirm button for downstream analysis.

Please note that the region must follow the pattern, `chromosome_id:start_site-end_site:strand`.

![](imgs/web/2.png)


At third tip, user could reset all configurations for another analysis. 
please be careful, this reset button will remove all information of previous plot.

### 3. Choose annotation

At first, user should click the `Refernce` option to select the genomic annotation file (GTF). And GTF without sorting or bgzipping are both fine for the tool.

Then user could define the parameter at the second tip, and click the confirm button to save the current information for the next step.
For choosing the annotation file, user could paste the absolute path of GTF or choose the file on webpage.
![](imgs/web/3.png)

Detailed explanation of configuration
![Detailed explanation of configuration](imgs/web/3.1.png)

Saved configurations
![](imgs/web/3.2.png)

### 4. Add plot by set input file

At first, user could choose the different type of plot for each dataset, and then prepare the parameter of each tract. 

![](imgs/web/4.png)

Please note that before processing another track, user should click the confirm button to save the current tract information.

![](imgs/web/4.1.png)

For another track, the document of parameter refer to [here](https://trackplot.readthedocs.io/en/latest/interactive/#api-documentation). 

### 5. Draw the final plot

After completing our configuration,user could define the parameter for output.

In addition to output the plot, the server provide preview option for user to quickly check the current configuration.

Once unexpected error happened, users could download the log file and create and issue on GitHub for us to debug.

![](imgs/web/5.png)

### 6. Preview image will display under the draw section

![](imgs/web/6.png)
