#!/usr/bin/env python

import os
import sys
from distutils.core import setup

sys.path.append(os.path.dirname(__file__)+"/python")
print sys.path

setup(name='ceo',
      version='1.0',
      description='Cuda--Engined Optics',
      author='Rodolphe Conan',
      author_email='conan.rod@gmail.com',
      url='http://rconan.github.io/CEO/',
      packages=['python.ceo'],
      install_requires=[
          'numpy','scipy',
      ])
