# -*- coding: utf-8 -*-

__version__ = "0.0.0"

try:
    __REBOUND_SETUP__
except NameError:
    __REBOUND_SETUP__ = False

if not __REBOUND_SETUP__:
    __all__ = []
