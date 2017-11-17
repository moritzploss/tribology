[![Build Status](https://travis-ci.org/moritzploss/tribology.png)](https://travis-ci.org/moritzploss/tribology)

# tribology
This package is a collection of methods and classes for tribology
research and education, including contact mechanics and lubrication
science. It provides implementations of analytical and numerical
calculation routines together with frequently used constants.

The package [documenation](https://moritzploss.github.io/tribology) is
under construction at the moment and provided through **GitHub Pages**
and in the `/docs` directory.

The **tribology** package is currently hosted on GitHub:

    https://github.com/moritzploss/tribology

You can also install the package through pip:

    pip install tribology

After installation, it is recommended to import the package and all its
modules as follows:

    ```python
    import tribology as tr
    ```

The package is developed by Moritz Ploss at KTH Royal
Institute of Technology, Stockholm, Sweden. For questions and comments,
please [email Moritz](mailto:moritz.ploss@gmail.com).

The package is provided under an MIT license. See the LICENSE.txt file
for more information.

# version history

#### 0.2.2.dev
Current development version on github
- docstrings for constants updated

#### 0.2.1
- Sphinx documentation added to git repo, docs available at https://moritzploss.github.io/tribology
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