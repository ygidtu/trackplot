#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""

"""


class AxLabel(object):

    __slots__ = ["Ax", "Label"]

    def __init__(self, ax, label):
        self.Ax = ax
        self.Label = label

    def __hash__(self):
        return hash(self.Ax)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


if __name__ == '__main__':
    pass
