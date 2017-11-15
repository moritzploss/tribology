from setuptools import setup

setup(
   name='tribology',
   version='0.1',
   description='a python 3 package for tribology research and education',
   author='Moritz Ploss',
   author_email='moritz.ploss@gmail.com',
   packages=['tribology'],
   install_requires=['numpy', 'numexpr', 'scipy'],
)
