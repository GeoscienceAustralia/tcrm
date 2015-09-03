The Tropical Cyclone Risk Model
===============================

The **Tropical Cyclone Risk Model** is a stochastic tropical cyclone 
model developed by
`Geoscience Australia <http://www.ga.gov.au>`_ for
estimating the wind hazard from tropical cyclones. 


Due to the relatively short record of quality-controlled, consistent tropical 
cyclone observations, it is difficult to estimate average recurrence interval 
wind speeds ue to tropical cyclones. To overcome the restriction of observed 
data, TCRM uses an autoregressive model to generate thousands of years of 
events that are statistically similar to the historical record. To translate 
these events to estimated wind speeds, TCRM applies a parametric windfield and 
boundary layer model to each event. Finally an extreme value distribution is 
fitted to the aggregated windfields at each grid point in the model domain to 
provide ARI wind speed estimates. 


Features
========


* **Multi-platform**: TCRM can run on desktop machines through to massively-parallel systems (tested on Windows XP/Vista/7, \*NIX);
* **Multiple options for wind field & boundary layer models**: A number of radial profiles and simple boundary layer models have been included to allow users to test sensitivity to these options.
* **Globally applicable**: Users can set up a domain in any TC basin in the globe. The model is not tuned to any one region of the globe. Rather, the model is designed to draw sufficient information from best-track archives;
* **Evaluation metrics**: Offers capability to run objective evaluation of track model metrics (e.g. landfall rates);
* **Single scenarios**: Users can run a single TC event (e.g. using a b-deck format track file) at high temporal resolution and extract time series data at chosen locations;



Dependencies
============

* TCRM requires `Python (2.7 preferred) <https://www.python.org/>`_,
  `Numpy <http://www.numpy.org/>`_, `Scipy <http://www.scipy.org/>`_,
  `Matplotlib <http://matplotlib.org/>`_, `Basemap
  <http://matplotlib.org/basemap/index.html>`_, `netcdf4-python
  <https://code.google.com/p/netcdf4-python/>`_, 
  `Shapely <https://github.com/Toblerity/Shapely>`_ and a C compiler;
* For parallel execution, `Pypar <http://github.com/daleroberts/pypar>`_ is required;

Status
======

.. image:: https://travis-ci.org/GeoscienceAustralia/tcrm.svg?branch=master
    :target: https://travis-ci.org/GeoscienceAustralia/tcrm
    :alt: Build status
    
.. image:: https://coveralls.io/repos/GeoscienceAustralia/tcrm/badge.svg?branch=master
  :target: https://coveralls.io/r/GeoscienceAustralia/tcrm?branch=master
  :alt: Test coverage

.. image:: https://landscape.io/github/GeoscienceAustralia/tcrm/master/landscape.svg?style=flat
   :target: https://landscape.io/github/GeoscienceAustralia/tcrm/master
   :alt: Code Health

Screenshot
==========

.. image:: https://rawgithub.com/GeoscienceAustralia/tcrm/master/docs/screenshot.png

Contributing to TCRM
====================

If you would like to take part in TCRM development, take a look at `docs/contributing.rst <https://github.com/GeoscienceAustralia/tcrm/blob/master/docs/contributing.rst>`_.

License information
===================

See the file `LICENSE.rst <https://github.com/GeoscienceAustralia/tcrm/blob/master/LICENSE.rst>`_ 
for information on the history of this software, terms and conditions for usage, 
and a DISCLAIMER OF ALL WARRANTIES.

