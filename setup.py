import os

from setuptools import setup, find_packages


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
    print(packages)
    return packages


setup(
    name='sashimi',
    long_description=__doc__,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    url='',
    entry_points={
        'console_scripts':
            [
                'sashimi = sashimi.cli:main'
            ]
    },
    python_requires='>=3.8',
    data_files=[(".", ['settings.ini'])],
    install_requires=locate_packages(),
)
