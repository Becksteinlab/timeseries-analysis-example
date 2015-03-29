====================================================================
Getting information out of molecular dynamics data: a simple example
====================================================================

This repository contains some example code illustrating how one
can write efficient analysis classes that collect timeseries data for each
frame of a molecular dynamics simulation trajectory. It is meant as
a pedagogical example, not as production code, and it is written to
illustrate the essential patterns to writing fast analyses.

The repository shows how to use this analysis class in an ipython notebook.
Note that this notebook requires at least ipython 3.0.0 to work. It also
requires pandas, pytables, and the hdf5 library. On Ubuntu 14.04 this means
using::

     apt-get install libhdf5-serial-1.8.4 libhdf5-serial-dev  

The notebook shows how to leverage using the efficient HDF5 file format with
pandas DataFrames, and gives a very short summary of some of the powerful
features of DataFrames for asking questions of your data.

