import os

from setuptools import setup, find_packages
from sashimi.cli import __version__


def locate_packages():
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    pipfile = os.path.join(__dir__, "Pipfile")
    packages = []
    if os.path.exists(pipfile):
        kept = False
        with open(pipfile) as r:
            for line in r:
                line = line.strip()

                if line == "[packages]":
                    kept = True
                    continue
                elif line.startswith("[") and line.endswith("]"):
                    kept = False

                if kept:
                    if line:
                        name, version = line.split(" = ")
                        version = version.strip('"')
                        if version == "*":
                            packages.append(f"{name}")
                        else:
                            packages.append(f"{name}, {version}")
    return packages


def load_description():
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    readme = os.path.join(__dir__, "README.md")
    if os.path.exists(readme):
        with open(readme) as r:
            return r.read()
    return "README"


setup(
    name='pysashimi',
    long_description=load_description(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    url="https://github.com/ygidtu/pysashimi",
    entry_points={
        'console_scripts':
            [
                'pysashimi = sashimi.cli:main'
            ]
    },
    python_requires='>=3.8',
    data_files=[(".", ['settings.ini'])],
    install_requires=locate_packages(),
    version=__version__
)
