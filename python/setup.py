#!/usr/bin/env python
# encoding: utf-8

import os
import sys
from functools import partial
from Cython.Build import cythonize

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension


# Publish the library to PyPI.
if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()


d = os.path.dirname
src_dir = os.path.join(d(d(os.path.abspath(__file__))), "src")

# Choose libraries to link.
libraries = []
if os.name == "posix":
    libraries.append("m")

# Specify the include directories.
include_dirs = [
    "rebound",
    src_dir,
]

rebound_sources = map(partial(os.path.join, src_dir), [
    "tree.c",
    "particle.c",
    "gravity_direct.c",
    "boundaries_open.c",
    "integrator_ias15.c",
    "input.c",
    "output.c",
    "collisions_none.c",
    "collision_resolve.c",
    "communication_mpi.c",
    "zpr.c",
    "display.c",
    "tools.c",
])

extensions = [
    Extension("rebound._rebound", libraries=libraries,
              include_dirs=include_dirs,
              sources=rebound_sources + [os.path.join("rebound",
                                                      "_rebound.pyx"),
                                         os.path.join("rebound",
                                                      "main.c")]),
]

# Hackishly inject a constant into builtins to enable importing of the
# package before the library is built.
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__REBOUND_SETUP__ = True
import rebound


# Execute the setup command.
desc = ""
# desc = open("README.rst").read()
setup(
    name="rebound",
    version=rebound.__version__,
    author="Daniel Foreman-Mackey",
    author_email="danfm@nyu.edu",
    packages=["rebound"],
    ext_modules=cythonize(extensions),
    # url="http://github.com/dfm/transit",
    # license="MIT",
    # description="A Python library for computing the light curves of "
    #             "transiting planets",
    # long_description=desc,
    # package_data={"": ["README.rst", "LICENSE", "include/*.h", ]},
    # include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
