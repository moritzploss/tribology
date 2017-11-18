.. Tribology documentation master file, created by
   sphinx-quickstart on Thu Nov 16 17:39:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tribology's Documentation!
=====================================

The tribology package is a collection of methods and classes for tribology
research and education, including contact mechanics and lubrication science.
The project is hosted on GitHub:

.. code-block:: python

   https://github.com/moritzploss/tribology

You can install the tribology package using pip:

.. code-block:: python

   pip install tribology

It is generally recommended to :code:`import` the tribology package as follows.
This will import all modules of the package and you can use all classes and
methods without ever worrying about modules:

.. code-block:: python

   import tribology as tr

Alternatively, you can perform module level imports for certain methods and
classes. Please see below for more information on available modules.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   sec_general_tribology.rst
   sec_contact_mechanics.rst
   sec_ehd_lubrication.rst
   sec_lubricants.rst
   sec_constants.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
