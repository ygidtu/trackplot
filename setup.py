import os
from setuptools import setup, find_packages

__version__ = "0.1.6"
__author__ = "ygidtu & Ran Zhou"
__email__ = "ygidtu@gmail.com"


def locate_packages():
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    pipfile = os.path.join(__dir__, "requirements.txt")

    packages = []
    with open(pipfile) as r:
        for line in r:
            line = line.strip()
            if not line.startswith("-") and \
                    "pybigwig" not in line.lower() and \
                    "hicmatrix" not in line.lower():
                packages.append(line.split(";")[0].strip())

    return packages


def load_description():
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    readme = os.path.join(__dir__, "README.md")
    if os.path.exists(readme):
        with open(readme) as r:
            return r.read()
    return "README"


setup(
    name='sashimi.py',
    long_description=load_description(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    url="https://github.com/ygidtu/sashimi.py",
    entry_points={
        'console_scripts': ['sashimipy=sashimi.cli:main']
    },
    python_requires='>=3.8',
    data_files=[(".", ['README.md', 'Pipfile', 'pyproject.toml', 'requirements.txt'])],
    install_requires=locate_packages(),
    extra_requires={
        "bigwig": ["pybigwig"],
        "hic": ["hicmatrix"]
    },
    version=__version__
)
