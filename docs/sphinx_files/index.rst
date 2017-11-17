.. Tribology documentation master file, created by
   sphinx-quickstart on Thu Nov 16 17:39:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tribology's Documentation!
=====================================

The documentation is currently under construction. The information below is
incomplete and not formated; this will be fixed soon. In the meantime, I hope
what you can find below is better than not having any documentation at all.

You can install the tribology package using pip:

.. code-block:: python

   pip install tribology

It is generally recommended to :code:`import` the tribology package as follows:

.. code-block:: python

   import tribology as tr

This will import all package modules. After that, you can use the below
methods and classes without the module name (numpy style). That means that the
two version of the following code have the same effect, but the first one is
probably more convenient:

.. code-block:: python

   # this code imports the beinflumat method from the
   # tribology_boundary_element module of the tribology package

   # without explicit module import
   import tribology as tr
   inf_mat = tr.beinflumat(...)

   # with explicit module import
   from tribology import tribology_boundary_element as tbe
   inf_mat = tbe.beinflumat(...)


.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   tribology.rst
   tribology_boundary_element.rst
   tribology_dowson_hamrock.rst
   tribology_hertz.rst
   tribology_lubrication.rst
   tribology_constants.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
