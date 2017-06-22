#!/usr/bin/env python

import os
import sys
import distutils.cmd
import distutils.log
import setuptools
import subprocess
from distutils.core import setup
import setuptools.command.build_py

sys.path.append(os.path.dirname(__file__)+"/python")
print sys.path

class MakeCeoCommand(distutils.cmd.Command):
  """A custom command to run Pylint on all Python source files."""

  description = 'Make CEO'
  user_options = [
      # The format is (long option, short option, description).
      ('none=', None, ''),
  ]

  def initialize_options(self):
    """Set default values for options."""
    # Each user option must be listed here with their default value.
    pass

  def finalize_options(self):
    """Post-process options."""
    pass

  def run(self):
    """Run command."""
    command = ['/usr/bin/make']
    #if self.pylint_rcfile:
    #  command.append('--rcfile=%s' % self.pylint_rcfile)
    #command.append(os.getcwd())
    command.append('all')
    command.append('cython')
    self.announce(
        'Running command: %s' % str(command),
        level=distutils.log.INFO)
    subprocess.check_call(command)

class BuildPyCommand(setuptools.command.build_py.build_py):
  """Custom build command."""

  def run(self):
    self.run_command('make_ceo')
    setuptools.command.build_py.build_py.run(self)


setuptools.setup(
    cmdclass={
        'make_ceo': MakeCeoCommand,
        'build_py': BuildPyCommand,
    },
    name='ceo',
    version='1.0',
    description='Cuda--Engined Optics',
    author='Rodolphe Conan',
    author_email='conan.rod@gmail.com',
    url='http://rconan.github.io/CEO/',
    packages=['python.ceo']
)
