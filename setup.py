from setuptools import setup

setup(
   name='tribology',
   version='0.5.13',
   setup_requires=['setuptools-git-version'],
   description='functions and constants for tribology research',
   long_description='A collection of functions and constants for tribology '
                    'research and education, including contact mechanics, '
                    'lubrication science, data processing and data analysis.',
   author='Moritz Ploss',
   author_email='moritz.ploss@gmail.com',
   packages=['tribology',
             'tribology.p3can',
             'tribology.tests'
             ],
   include_package_data=True,
   license='MIT',
   keywords='tribology machine design research',
   install_requires=[
       'numpy<=1.16.1',
       'numexpr',
       'scipy',
       'opencv-python==3.4.2.16',
       'matplotlib',
   ],
   python_requires='>=3.5',
   url='https://github.com/moritzploss/tribology',
   classifiers=[
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 4 - Beta',

    # Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',

    # Pick your license as you wish (should match "license" above)
    'License :: OSI Approved :: MIT License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6'],
)
