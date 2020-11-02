Examples
========

Full TCRM hazard simulation
---------------------------

An example configuration file is provided with the code to provide
users a starting point for testing the installation, and for setting
up their own simulation. The example simulation is centred on Port
Hedland, Australia (118.6E, 20.3S), which has experienced numerous
severe tropical cyclones.

Once the model has been :ref:`installed <installation>` and tested,
the example simulation can be run as follows::
    
    $ python tcrm.py -c example/port_hedland.ini

The model will automatically create output folders to store the
results in, as well as a log file that provides details on the
progress of the model. The simulation will process the input track
database (IBTrACS v03r10), run the statistical interface, generate
1000 simulated years of TC events, the wind swaths associated with
those simulated years, calculate the hazard (return period wind
speeds) for the region and plot the output as maps and a return period
curve for Port Hedland. The hazard data are also stored in netCDF
(version 4) files.

.. figure:: /docs/images/hazard_example.png
     :align: center
     :alt: 100-year return period wind speed near Port Hedland,
           Australia.
     :figclass: align-center

     100-year return period wind speed near Port Hedland,
     Australia. Wind speeds represent a 0.2-second gust wind speed, 10
     metres above ground level in open, flat terrain.

.. figure:: /docs/images/hazard_curve.png
    :align: center
    :alt: Example hazard curve for Port Hedland, Australia.
    :figclass: align-center
    
    Hazard curve for Port Hedland, Australia, based on 1000 years of
    synthetic tropical cyclone events. The red points are the wind speeds for
    individual simulated events, and the blue line is a fitted generalised
    pareto distribution (GPD).

Running in an MPI environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For a multiprocessor environment using mpi4py, the example can
be run with::

    $ mpirun -np 16 python tcrm.py -c example/port_hedland.ini

This will execute TCRM across 16 CPUs (if available). 

TCRM has been tested on systems up to 256 CPUs. Testing with moderate
event sets (4000 events) indicate > 23 times speedup when run across
32 CPUs.

Scenario simulation
-------------------

An example scenario is also included with the code to demonstrate an
individual event simulation. This uses Tropical Cyclone *Yas*, which
impacted the north Quensland coast in February 2011 as a severe
category 4 cyclone. The example simulation uses best-track data from
the Bureau of Meteorology. See the :ref:`scenariomodelling` section
for more details.


