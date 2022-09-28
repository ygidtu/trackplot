#!/usr/bin/env python3
# -*-coding:utf-8 -*-

u"""
Handle the customized junction
"""
import os
import re
from typing import Dict

from sashimi.base.Junction import Junction


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
                    match = re.search(r"(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+)", line[0])

                    if match:
                        match = match.groupdict()

                        junc = Junction(match["chrom"], int(match["start"]), int(match["end"]))

                        for i in range(1, len(header)):
                            sample = header[i]
                            if sample not in data.keys():
                                data[sample] = {}
                            data[sample][junc] = int(line[i])

    return data


if __name__ == '__main__':
    pass
