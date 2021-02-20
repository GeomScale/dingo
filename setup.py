
# This is the setup Python script for building the volestipy Python module

from distutils.core import setup
from distutils.core import Extension
from Cython.Build import cythonize
from os.path import join
import numpy
import os
import git

# TODO: update these information
version = "0.3.0"
license='LGPL3',
packages = ["volestipy"]
description="volestipy: wrapper for the VolEsti library to sample from convex sets and compute volume."
author = "Pedro Zuidberg Dos Martires, Haris Zafeiropoulos"
author_email="pedro.zuidbergdosmartires@cs.kuleuven.be, haris-zaf@hcmr.gr"
name = 'volestipy'
#zip_safe = False

if (not os.path.isdir("eigen")) :
    print('Cloning eigen library...')
    git.Repo.clone_from('https://gitlab.com/libeigen/eigen.git', 'eigen', branch='3.3')

source_directory_list = ['volestipy', join('volestipy','src')]

compiler_args = [
 "-std=c++11",
 "-O3",
  "-DBOOST_NO_AUTO_PTR",
 "-ldl",
 "-lm"
]

link_args = ['-O3']

extra_volesti_include_dirs = [
# inside the volestipy directory, there is an extra volestipy folder containing
# a "include" folder; we add those as included dirs; the bindings.h file is located there
 join("volestipy","include"),

# the volesti code uses some external classes. these are located on the "external"
# directory and we need to add them as well
# join("..","additional_external"),
 join("eigen"),
 join("volesti","external"),
 join("volesti","external","minimum_ellipsoid"),
# join("..","volesti","external","LPsolve_src","run_headers"),
 join("volesti","external","boost"),

# we also move back and include and add the directories on the "include" directory
# (generatorors, random_walks, sampling etc)
 join("volesti","include"),
 join("volesti","include","convex"),
 join("volesti","include","misc"),
 join("volesti","include","random_walks"),
 join("volesti","include","volume"),
 join("volesti","include","generators"),
 join("volesti","include","cartesian_geom"),
]

src_files = ["volestipy/volestipy.pyx","volestipy/src/bindings.cpp"]
extra_include_dirs = [numpy.get_include()]
# Return the directory that contains the NumPy *.h header files.
# Extension modules that need to compile against NumPy should use this function to locate the appropriate include directory.

# here I need to check how it is actually included this library
#library_includes = ["lpsolve55"]

ext_module = Extension(
 "volestipy",
 language = "c++",
 sources = src_files,
 include_dirs = extra_include_dirs + extra_volesti_include_dirs,
 extra_compile_args = compiler_args,
 extra_link_args = link_args,
)
print("The Extension function is OK.")

ext_modules = cythonize([ext_module], gdb_debug=False)
print("The cythonize function ran fine!")

setup(
 version = version,
 author = author,
 author_email = author_email,
 name = name,
 packages = packages,
 ext_modules = ext_modules
# zip_safe=zip_safe,
)

print("Finally, the setup function was performed. We ve got out interface.")
