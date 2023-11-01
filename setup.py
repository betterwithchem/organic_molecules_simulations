from setuptools import setup

setup(
   name='sim_launch_py',
   version='0.1a',
   description='A module to prepare and run MD simulations',
   author='Matteo Paloni',
   author_email='m.paloni@ucl.ac.uk',
   packages=['sim_launch_py'],  #same as name
   install_requires=['numpy'], #external packages as dependencies
)
