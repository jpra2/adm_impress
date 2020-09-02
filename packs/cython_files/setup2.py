from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = Extension(name='test_ent', sources=['test_ent.pyx'])
setup(
    ext_modules=cythonize(extensions, language_level=3)
)

# python setup.py build_ext --inplace
