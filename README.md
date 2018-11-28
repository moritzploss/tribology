[![Build Status](https://travis-ci.org/moritzploss/tribology.png)](https://travis-ci.org/moritzploss/tribology)
[![DOI](https://zenodo.org/badge/110825481.svg)](https://zenodo.org/badge/latestdoi/110825481)


# tribology
This Python 3 package is a collection of functions for
tribology research and education, including contact mechanics,
lubrication science, data handling and data processing. It provides 
implementations of analytical and numerical calculation routines together with 
frequently used constants.

The **tribology package** is hosted on GitHub:

    https://github.com/moritzploss/tribology

You can **install** the package from the
**<a href="https://pypi.python.org/pypi/tribology" target="_blank">PyPI</a>**
index using pip:

    pip install tribology

After installation, it is recommended to **import** the package and all
its modules as follows:

```python
import tribology as tr
```

The package
**<a href="https://moritzploss.github.io/tribology" target="_blank">
documenation</a>** is provided through GitHub Pages; the Sphinx source files
can be found in the [`/docs`](./docs) directory. Simple examples of how to use 
the package are provided in the [`/examples`](./examples) directory.

The package is provided under an MIT license. See the 
[`LICENSE.txt`](LICENSE.txt) file for more information.


# known issues
Some older versions of the package (0.4 <= version < 0.5) may raise the 
following error when imported:

    ModuleNotFoundError: No module named 'cv2'

This is caused by a missing `install_requires` argument in
[`setup.py`](./setup.py). To fix this, either update to a later version of the
package, or install the missing `cv2` package manually:

    pip install opencv-python


# use in scientific publications

You can refer to the tribology package in **scientific publications** by
using its DOI. The following DOI will resolve all releases of the
package and automatically point to the latest release. A more detailed
overview of the DOI release history can be found through
**<a href="https://doi.org/10.5281/zenodo.1117727" target="_blank">Zenodo</a>**:

    DOI:  10.5281/zenodo.1117727

The suggested citation format is:

    [List of contributors] et al., Tribology -- A Python Package for Tribology 
    Research and Education, version v[x.x.x] ([date]). doi: [DOI] 


# version history

Please refer to the release history for a detailed overview. The current version
on GitHub incorporates the following changes with respect to the latest release:

- None
