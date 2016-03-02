#!/usr/bin/env python

import sys
from distutils.core import setup

sys.path.append("python")

setup(name='ceo',
      version='1.0',
      description='Cuda--Engined Optics',
      author='Rodolphe Conan',
      author_email='conan.rod@gmail.com',
      url='http://rconan.github.io/CEO/',
      packages=['python.ceo'],
      install_requires=[
          'numpy',
      ])
