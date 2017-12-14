[![Build Status](https://travis-ci.org/moritzploss/tribology.png)](https://travis-ci.org/moritzploss/tribology)

# tribology
This Python 3 package is a collection of methods and classes for
tribology research and education, including contact mechanics,
lubrication science and data handling. It provides implementations of
analytical and numerical calculation routines together with frequently
used constants.

The **tribology** package is currently hosted on GitHub:

    https://github.com/moritzploss/tribology

You can also install the package through pip:

    pip install tribology

After installation, it is recommended to import the package and all its
modules as follows:

```python
import tribology as tr
```

The package
**<a href="https://moritzploss.github.io/tribology" target="_blank">
documenation</a>** is provided through GitHub Pages and can be found
in the `/docs` directory. Simple examples of how to use the package are
provided in the `/demos`
directory.

The package is developed by Moritz Ploss at KTH Royal
Institute of Technology, Stockholm, Sweden. For questions and comments,
please **[send an email to Moritz](mailto:moritz.ploss@gmail.com)**.

The package is provided under an MIT license. See the `LICENSE.txt` file
for more information.

# version history

Neutral builds are continuously deployed to PyPI. The log for
neutral builds includes changes with respect to the latest numbered
release listed below.

#### Latest Neutral Build

- method `walther` now supports `ndarray` arguments
- method `profrolleriso` added to `tribology` module
- no longer compatible with Python 3.4 and lower
- the `data_import` module now provides methods to import output files
from PCS Instruments test rigs.
- module `data_import` added. the modgit gitule contains methods for data
import from delimited data files into Numpy and Matlab database format.
- modules renamed, removed leading `tribology_` in module names
- auto-release on PyPI through Travis CI
- refactoring

#### 0.2.2

- html documentation now with multi-page structure
- docstrings changed to numpy format and refactored
- functions `eeff` and `reff` moved to `tribology_hertz` module
- function `eeff` renamed to `reff`
- function `meff` renamed to `eeff`
- helper functions made private

#### 0.2.1
- Sphinx documentation added to GitHub repo, docs available at
https://moritzploss.github.io/tribology
-  Travis CI builds now with Pylint error check for package files.
Non-package files are not checked. Build fails if Pylint error is found.
- Python 3.4 and 3.5 builds added to Travis CI (now 3.4, 3.5 and 3.6)
- change log added to `README.md`
- various methods renamed
- helper functions made private
- classifiers added to `setup.py`
- demo for boundary element methods added (`demo_ball_on_plate.py`)
- refactoring

#### 0.2.0
first release on PyPI