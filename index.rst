.. TCRM documentation master file, created by
   sphinx-quickstart on Mon Apr 07 16:22:26 2014.
   You can adapt this file completely to your liking,
   but it should at least contain the root `toctree` 
   directive.

The Tropical Cyclone Risk Model
===============================

The **Tropical Cyclone Risk Model** is a stochastic tropical cyclone 
model developed by
`Geoscience Australia <http://www.ga.gov.au>`_ for
estimating the wind hazard from tropical cyclones. 

Due to the relatively short record of quality-controlled, consistent
tropical cyclone observations, it is difficult to estimate average
recurrence interval wind speeds due to tropical cyclones. To overcome
the restriction of observed data, TCRM uses an autoregressive model to
generate thousands of years of events that are statistically similar
to the historical record. To translate these events to estimated wind
speeds, TCRM applies a parametric windfield and boundary layer model
to each event. Finally an extreme value distribution is fitted to the
aggregated windfields at each grid point in the model domain to
provide ARI wind speed estimates.

Features
========

* **Multi-platform**: TCRM can run on desktop machines through to 
  massively-parallel systems (tested on Windows XP/Vista/7, \*NIX);
* **Multiple options for wind field & boundary layer models**: A 
  number of radial profiles and simple boundary layer models have 
  been included to allow users to test sensitivity to these options.
* **Globally applicable**: Users can set up a domain in any TC basin 
  in the globe. The model is not tuned to any one region of the 
  globe. Rather, the model is designed to draw sufficient 
  information from best-track archives or TC databases;
* **Evaluation metrics**: Offers capability to run objective 
  evaluation of track model metrics (e.g. landfall rates);
* **Single scenarios**: Users can run a single TC event (e.g. using 
  a b-deck format track file) at high temporal resolution and 
  extract time series data at chosen locations;

Releases
========

Latest releases can be downloaded from the `Geoscience Australia GitHub repository <https://github.com/GeoscienceAustralia/tcrm/releases>`_. 

Bleeding edge versions are accessible `here <https://github.com/GeoscienceAustralia/tcrm/archive/master.zip>`_. 

Contributions are welcome -- create a `fork
<https://github.com/GeoscienceAustralia/tcrm/fork>`_ or clone the
repo.

Contents
========
.. toctree::
   :maxdepth: 1

   Introduction <docs/intro>
   Installation <docs/install>
   Setting up the model <docs/setup>
   Running the model <docs/execution>
   Examples <docs/examples>
   Model evaluation <docs/evaluation>
   Using different TC databases <docs/sources>
   Scenario simulation <docs/scenarios>
   Contributing to TCRM <docs/contributing>
   Utilities <docs/Utilities>
   Glossary <docs/glossary>
   License information <LICENSE>

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

