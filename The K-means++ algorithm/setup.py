from setuptools import setup, Extension, find_packages

setup(name='mykmeanssp',
      version='1.0',
      description='HW2 maya ali',
      ext_modules=[Extension('mykmeanssp', sources=['kmeansmodule.c'])])



