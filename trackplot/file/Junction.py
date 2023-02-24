#!/usr/bin/env python3
# -*-coding:utf-8 -*-

u"""
Handle the customized junction
"""
import os
import re
from typing import Dict

from trackplot.base.Junction import Junction


def load_custom_junction(path: str) -> Dict[str, Dict[Junction, int]]:
    data = {}
    if os.path.exists(path) and os.path.isfile(path):
        header = []
        with open(path) as r:
            for line in r:
                if line.startswith("#"):
                    continue

                line = line.strip().split()

                if not header:
                    header = line
                else:
                    junc = Junction.create_junction(line[0])
                    for i in range(1, len(header)):
                        sample = header[i]
                        if sample not in data.keys():
                            data[sample] = {}
                        data[sample][junc] = int(line[i])

    return data


if __name__ == '__main__':
    pass
