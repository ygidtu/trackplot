#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
Backward-compatibility shim. Previously this file contained all plot rendering functions.
Now everything is in the ``trackplot.plot`` package (see trackplot/plot/).

This module re-exports all public symbols so that existing imports continue to work.
"""

from trackplot.plot import *  # noqa: F403
