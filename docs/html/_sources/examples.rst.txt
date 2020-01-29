
Examples
========

Compute a Single Transport Coefficient
--------------------------------------
To compute a single transport coefficient do the following

.. code-block:: python

   Am = 1.9944235e-23 # atomic mass of carbon-12 [g]
   rho_i = 1 # mass density [g/cc]
   T = 0.2 # temeprature [eV]
   Z = 6 # atomic number of carbon

   # Instantiate the stanton-murillo transport class
   sm = sm_transport(Am, rho_i, T, Z, units_out='cgs')

   # Compute transport coefficients
   D = sm.self_diffusion()
   eta = sm.viscocity()
   K = sm.thermal_conductivity()