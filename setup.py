from setuptools import setup

setup(
   name='tribology',
   version='0.1.24',
   description='methods and classes for tribology research',
   long_description='collection of methods and classes for tribology research '
                    'and education, including contact mechanics and '
                    'lubrication science,',
   author='Moritz Ploss',
   author_email='moritz.ploss@gmail.com',
   packages=['tribology'],
   license='MIT',
   keywords='tribology machine design research',
   install_requires=['numpy', 'numexpr', 'scipy'],
   python_requires='>=3.5',
   url='https://github.com/moritzploss/tribology'
)
