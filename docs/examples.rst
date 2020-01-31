.. role:: python(code)
   :language: python

Examples
========

Stanton and Murillo Transport Class
-----------------------------------
The Stanton and Murillo transport class (:python:`sm_transport()`) class allows users to compute single 
a transport coefficient or across a range of elements, mass densities, and temperatures using the 
`Stanton and Murillo model <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.93.043203>`_. 
Refer to *Fig 1.* for the structure of the output in the multi-element/mass-density/temperature case.

Computing a Single Transport Coefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To compute a single transport coefficient do the following:

.. code-block:: python

   import plasma_properties as pprops

   Am = 1.9944235e-23 # Atomic mass of element [g]
   rho_i = 1 # Mass density [g/cc]
   T = 0.2 # Temperature [eV]
   Z = 6 # Atomic number

   # Instantiate the Stanton-Murillo transport class
   sm = pprops.sm_transport(Am, rho_i, T, Z, units_out='cgs')
   
   # Compute transport coefficients
   D = sm.self_diffusion()
   eta = sm.viscocity()
   K = sm.thermal_conductivity()

   print(D, eta, K)

Output (in cgs units):

:code:`[0.00127858] [0.00084193] [21856.10137931]`


Computing Multiple Transport Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Computing transport coefficients for a range of elements, mass densities, and temperatures follows
the same syntax as the single transport case except the inputs are now :python:`numpy.array()`.

Below is example input for this case:

.. code-block:: python

   import plasma_properties as pprops

   Am = np.array([1.9944235e-23, 4.4803895e-23, 8.4590343e-23]) # Atomic masses for each element [g]
   rho_i = np.array([1,10,100]) # Mass densities [g/cc]
   T = np.arange(0.2, 200, 0.1) # Temperature range [eV]
   Z = np.array([6, 13, 23]) # Atomic number for each element

   # Instantiate the Stanton-Murillo transport class
   sm = sm_transport(Am, rho_i, T, Z, units_out='cgs')

   # Compute transport 
   D = sm.self_diffusion()
   eta = sm.viscocity()
   K = sm.thermal_conductivity()

The output of the above code will be a 3-dimensional :python:`numpy.ndarray()` where the axes have the following structure.

.. figure:: _images/transport_data_structure_grid2.png
 :width: 400
 :align: center
 :alt: Transport Coefficient Data Structure

 Fig 1. The shape of the data structure that is output for the case of multi-element/temperature/density transport coefficients. 
 Note that each 2-dimensional "slice" in the *Z* direction corresponds to a different element, and moving along the positive :math:`\rho`/T direction corresponds to an increase in the mass-density/temperature for a fixed element.

.. note::

   Referencing this data structure directly corresponds to the entries in each of the input arrays. For example, if you wish to print the self-diffusion coefficient from the above code block for carbon (0th element of the *Z* array), at 10 g/cc (1st element of the *rho_i* array), at 0.4 eV (2nd element of the *T* array), you would use the syntax
   :python:`print(D[1,2,0])` (marked in red in *Fig. 1*). 

Below is some example code for plotting data in the multiple element/mass-density/temperature case:

.. code-block:: python
   
   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(1, 3, figsize=(30,8))

   #---------------- Plotting Self-Diffusion ----------------#
   ax[0].loglog(T, D[0,:,0], linewidth=3, label='Carbon')
   ax[0].loglog(T, D[0,:,1], linewidth=3, label='Aluminum')
   ax[0].loglog(T, D[0,:,2], linewidth=3, label='Vanadium')

   ax[0].set_xlabel('Temperature [eV]', fontsize=20)
   ax[0].set_ylabel('Self-Diffusion $[cm^2/s]$', fontsize=20)
   ax[0].legend(fontsize=18)
   ax[0].tick_params(axis="x", labelsize=18) 
   ax[0].tick_params(axis="y", labelsize=18) 


   #------------------ Plotting Viscosity -------------------#
   ax[1].loglog(T, eta[0,:,0], linewidth=3, label='Carbon')
   ax[1].loglog(T, eta[0,:,1], linewidth=3, label='Aluminum')
   ax[1].loglog(T, eta[0,:,2], linewidth=3, label='Vanadium')

   ax[1].set_xlabel('Temperature [eV]', fontsize=20)
   ax[1].set_ylabel('Viscosity $[g/(cm * s)]$', fontsize=20)
   ax[1].legend(fontsize=18)
   ax[1].tick_params(axis="x", labelsize=18) 
   ax[1].tick_params(axis="y", labelsize=18) 


   #-------------- Plotting Thermal Conductivity ------------#
   ax[2].loglog(T, K[0,:,0], linewidth=3, label='Carbon')
   ax[2].loglog(T, K[0,:,1], linewidth=3, label='Aluminum')
   ax[2].loglog(T, K[0,:,2], linewidth=3, label='Vanadium')

   ax[2].set_xlabel('Temperature [eV]', fontsize=20)
   ax[2].set_ylabel('Thermal Conductivity $[erg/(cm * s * K)]$', fontsize=20)
   ax[2].legend(fontsize=18)
   ax[2].tick_params(axis="x", labelsize=18) 
   ax[2].tick_params(axis="y", labelsize=18) 

   plt.show()
   plt.savefig('transport_compare.png', dpi=300, bbox_inches='tight')

.. figure:: _images/transport_compare.png
 :width: 850
 :align: center
 :alt: Self-diffusion, viscosity, and thermal conductivity plots


Example: Viscosity versus Density
*********************************

.. code-block:: python

   import plasma_properties as pprops
   import matplotlib.pyplot as plt

   Am = np.array([1.9944235e-23, 4.4803895e-23, 8.4590343e-23]) # atomic masses for each element [g]
   rho_i = np.arange(0.1, 10, 0.01) # mass density [g/cc]
   T = 0.2 # temeprature [eV]
   Z = np.array([6, 13, 23]) # atomic number for each element

   # Instantiate the stanton-murillo transport class
   sm = pprops.sm_transport(Am, rho_i, T, Z, units_out='cgs')

   # Compute viscocity
   eta = sm.viscocity()

   # Plotting
   plt.figure(figsize=(10,8))

   # Plot viscosity versus density for fixed temp (0.1 eV)
   plt.loglog(rho_i, eta[:,0,0], linewidth=3, label='Carbon')
   plt.loglog(rho_i, eta[:,0,1], linewidth=3, label='Aluminum')
   plt.loglog(rho_i, eta[:,0,2], linewidth=3, label='Vanadium')

   plt.xlabel('Mass-Density [g/cc]', fontsize=20)
   plt.ylabel('Viscosity $[g/(cm s)]$', fontsize=20)
   plt.legend(fontsize=18, loc='upper left')
   plt.yticks([7e-4, 1e-3, 2e-3])
   plt.tick_params(axis="x", labelsize=16) 
   plt.tick_params(axis="y", labelsize=16)
   plt.show()
   plt.savefig('viscosity.png', dpi=300, bbox_inches='tight')

.. figure:: _images/viscosity.png
 :width: 400
 :align: center
 :alt: viscosity coefficient as a function of density
