#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
This file contains the configuration of different themes
"""
from matplotlib import axes


class Theme(object):
    @classmethod
    def blank(cls, ax: axes.Axes):
        u"""

        :param ax:
        :return:
        """
        # ax.spines['left'].set_position('zero')
        # ax.spines['bottom'].set_position('zero')
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        ax.set_axis_off()
        ax.tick_params(bottom=False, top=False, left=False, right=False)

    @classmethod
    def ticks(cls, ax: axes.Axes):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(bottom=True, top=False, left=True, right=False)

    @classmethod
    def ticks_blank(cls, ax: axes.Axes):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticklabels([])
        ax.tick_params(bottom=False, top=False, left=True, right=False)

    @classmethod
    def set_theme(cls, ax: axes.Axes, name: str = "blank"):
        for i, func in cls.__dict__.items():
            if isinstance(func, classmethod) and i != "get" and i == name:
                """
                This usage see: https://blog.peterlamut.com/tag/classmethod/
                """
                return func.__get__(None, cls)(ax)

        raise ValueError(f"{name} is not a valid theme")


if __name__ == "__main__":
    pass
