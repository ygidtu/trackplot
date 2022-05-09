#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""

"""
from typing import List


class Stroke(object):
    def __init__(self, start: int, end: int, color: str = "red", label: str = ""):
        self.start = start
        self.end = end
        self.color = color
        self.label = label

    @property
    def center(self) -> float:
        return (self.end + self.start) / 2

    @classmethod
    def create_from_string(cls, stroke: str):
        res = []
        for i in stroke.split(":"):
            i = i.split("@")
            sites = sorted([int(x) - self.start for x in i[0].split("-")])
            if sites[0] < 0:
                sites[0] = 0

            if sites[-1] > len(self):
                sites[-1] = len(self)

            color = "red"
            label = ""
            if len(i) > 1:
                i = i[-1].split("-")
                color = i[0]

                if len(i) > 1:
                    label = i[-1]

            res.append(Stroke(sites[0], sites[-1], color, label))
        return res


if __name__ == "__main__":
    pass
