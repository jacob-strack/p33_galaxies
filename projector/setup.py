#setup script for cython loop 

from distutils.core import setup
from Cython.Build import cythonize

setup(
        ext_modules=cythonize("projector/loop.pyx", compiler_directives={"language_level":"3"}
            )
    )
