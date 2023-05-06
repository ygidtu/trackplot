import os
from configparser import ConfigParser
from setuptools import setup, find_packages

__version__ = "0.2.4"
__author__ = "ygidtu & Ran Zhou"
__email__ = "ygidtu@gmail.com"


def locate_packages():
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    pipfile = os.path.join(__dir__, "Pipfile")

    conf = ConfigParser()
    conf.read(pipfile)

    packages = []
    for name in conf["packages"]:
        version = conf.get("packages", name)
        version = version.strip('"').strip("'").lstrip("^")
        if version == "*":
            packages.append(name)
        else:
            packages.append(f"{name}=={version}")
    return packages


def load_description():
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    readme = os.path.join(__dir__, "README.md")
    if os.path.exists(readme):
        with open(readme) as r:
            return r.read()
    return "README"


setup(
    name='trackplot',
    long_description=load_description(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    url="https://github.com/ygidtu/trackplot",
    entry_points={
        'console_scripts':
            [
                'trackplot = trackplot.cli:main'
            ]
    },
    python_requires='>=3.8',
    data_files=[(".", ['README.md', 'Pipfile', 'Pipfile.lock'])],
    install_requires=locate_packages(),
    version=__version__
)
