import os

# See if Cython is installed
try:
    from Cython.Build import cythonize
# Do nothing if Cython is not available
except ImportError:
    # Got to provide this function. Otherwise, poetry will fail
    def build(setup_kwargs):
        pass


# Cython is installed. Compile
else:
    from setuptools import Extension
    from setuptools.dist import Distribution
    from distutils.command.build_ext import build_ext

    # This function will be executed in setup.py:
    def build(setup_kwargs):
        # The file you want to compile
        extensions = ["dingo/volestipy.pyx"]

        # gcc arguments hack: enable optimizations
        os.environ["CFLAGS"] = [
            "-std=c++11",
            "-O3",
            "-DBOOST_NO_AUTO_PTR",
            "-ldl",
            "-lm",
        ]

        # Build
        setup_kwargs.update(
            {
                "ext_modules": cythonize(
                    extensions,
                    language_level=3,
                    compiler_directives={"linetrace": True},
                ),
                "cmdclass": {"build_ext": build_ext},
            }
        )
