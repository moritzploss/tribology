.. Tribology documentation master file, created by
   sphinx-quickstart on Thu Nov 16 17:39:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tribology's Documentation!
=====================================

This Python 3 package is a collection of functions and classes for tribology
research and education, including contact mechanics, lubrication
science and data handling. The project is hosted on GitHub:

.. code-block:: python

   https://github.com/moritzploss/tribology

You can install the tribology package using pip:

.. code-block:: python

   pip install tribology

You can refer to the tribology package in **scientific publications** by
using its DOI. The following DOI will resolve all releases of the
package and automatically point to the latest release:

.. code-block:: python

   DOI: 10.5281/zenodo.1117727

It is generally recommended to :code:`import` the tribology package as follows.
This will import all modules of the package and you can use all classes and
functions without ever worrying about modules:

.. code-block:: python

   import tribology as tr

Alternatively, you can perform module level imports using the moule names given
in the remainder of this documentation. For quick access to relevant
information, package modules are loosely grouped by topic in the following.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   sec_general_tribology.rst
   sec_contact_mechanics.rst
   sec_ehd_lubrication.rst
   sec_lubricants.rst
   sec_data_import.rst
   sec_constants.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
