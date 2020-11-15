.. _single-transport:


========================================
Computing A Single Transport Coefficient
========================================

This is an example of using `~plasma_properties.transport.SM~.

.. code-block:: python

   from plasma_properties import transport

   Am = 1.9944235e-23 # Atomic mass of element [g]
   rho_i = 1 # Mass density [g/cc]
   T = 0.2 # Temperature [eV]
   Z = 6 # Nuclear charge for carbon

   # Instantiate the Stanton-Murillo transport submodule
   sm = transport.SM(Am, rho_i, T, Z, units_out='cgs')
   
   # Compute transport coefficients
   D = sm.self_diffusion()
   eta = sm.viscosity()
   K = sm.thermal_conductivity()

   print(D, eta, K)

Output (in cgs units):

:code:`[0.00127858] [0.00084193] [21856.10137931]`