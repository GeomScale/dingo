# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

# This is the setup Python script for building the dingo library

from distutils.core import setup
from distutils.core import Extension
from Cython.Build import cythonize
from os.path import join
import numpy
import os

# information about the dingo library
version = "0.1.0"
license = ("LGPL3",)
packages = ["dingo"]
description = "A python library for metabolic networks sampling and analysis"
author = "Apostolos Chalkis"
author_email = "tolis.chal@gmail.com"
name = "dingo"


source_directory_list = ["dingo", join("dingo", "bindings")]

compiler_args = ["-std=c++11", "-O3", "-DBOOST_NO_AUTO_PTR", "-ldl", "-lm", "-fopenmp"]

link_args = ["-O3", "-fopenmp"]

extra_volesti_include_dirs = [
    # include binding files
    join("dingo", "bindings"),
    # the volesti code uses some external classes.
    # external directories we need to add
    join("eigen"),
    join("boost_1_76_0"),
    join("boost_1_76_0", "boost"),
    join("volesti", "external"),
    join("volesti", "external", "minimum_ellipsoid"),
    # include and add the directories on the "include" directory
    join("volesti", "include"),
    join("volesti", "include", "convex_bodies"),
    join("volesti", "include", "random_walks"),
    join("volesti", "include", "volume"),
    join("volesti", "include", "generators"),
    join("volesti", "include", "cartesian_geom"),
]

src_files = ["dingo/volestipy.pyx", "dingo/bindings/bindings.cpp"]

# Return the directory that contains the NumPy *.h header files.
# Extension modules that need to compile against NumPy should use this
# function to locate the appropriate include directory.
extra_include_dirs = [numpy.get_include()]

ext_module = Extension(
    "volestipy",
    language="c++",
    sources=src_files,
    include_dirs=extra_include_dirs + extra_volesti_include_dirs,
    extra_compile_args=compiler_args,
    extra_link_args=link_args,
)
print("The Extension function is OK.")

ext_modules = cythonize([ext_module], gdb_debug=False)
print("The cythonize function ran fine!")

setup(
    version=version,
    author=author,
    author_email=author_email,
    name=name,
    packages=packages,
    ext_modules=ext_modules,
)

print("Installation of dingo completed.")
