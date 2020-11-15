import numpy as np
import matplotlib.pyplot as plt
from plasma_properties import zbar
from plasma_properties import error
import parameters

class SM:
    """Generate the Stanton and Murillo transport coefficients. test
    
    Parameters
    ----------
    Am : float or arrary_like
        Atomic mass of element (or isotope) in units of grams [g].

    mass_density : float or array_like
        Range of mass densities in units of grams per 
        cubic centimeter [g/cc]. 
    T : float or array_like
        Temperature range in units of electron-volts [eV]

    Z : int or arrray_like
        Atomic number for each element
   
    units_out : str
        Unit system for resulting transport coefficient.
        Default: dimensionless "star" units.

    References
    ----------
    .. [1] `Stanton, Liam G., and Michael S. Murillo. "Ionic transport in high-energy-density matter." Physical Review E 93.4 (2016): 043203.
        <https://doi.org/10.1103/PhysRevE.93.043203>`_
    """

    def __init__(self, Am, mass_density, T, Z, units_out='star'):
        """ Initialize all parameters to be shared across methods.
        """
        for X in [Am, mass_density, T, Z]:
            if str(type(X)) == "<class 'list'>":
                raise TypeError('type "list" not supported')

        # Check type of input and deal with float cases
        self.mi = Am
        if str(type(self.mi)) != "<class 'numpy.ndarray'>":
            self.mi = np.array([Am])

        self.Z = Z
        if str(type(self.Z)) != "<class 'numpy.ndarray'>":
            self.Z = np.array([Z])
       
        if str(type(mass_density)) != "<class 'numpy.ndarray'>":
            self.num_density = np.array([mass_density])/Am
            print
        else:
            if len(mass_density) > 1:
                self.num_density = np.tile(mass_density, len(Am))/np.repeat(Am, len(mass_density))
                self.num_density = np.reshape(self.num_density, (len(Am) ,len(mass_density)))
            else:
                self.num_density = mass_density/Am
                
        self.T = T
        if str(type(self.T)) != "<class 'numpy.ndarray'>":
            self.T = np.array([T])

        # Check input, throw errors
        if len(self.mi) != len(self.Z):
            raise DimensionMismatchError('atomic mass and nuclear charge arrays must be the same size')

        # for X, st in [(self.mi,'Am'), (self.num_density, 'mass_density'), (self.T, 'T'), (self.Z, 'Z')]:
        #     print(X)
        #     if any(v <= 0 for v in X) :
        #         raise ValueError('{} values must be > 0'.format(st))

        # Class wide parameter definition in cgs units
        self.e_squared = 1.4399644800e-7 # [ev cm]
        self.hbar = 6.5822958e-16 # reduced Planck's constant [eV*s]
        self.me = 9.1095e-28 # mass of electron [g]
        self.units_out = units_out
        print('Transport coefficients returned in {} units.'.format(self.units_out))

    def plasma_params(self, Z, ni, mi):
        """ Computes the plasma parameters (e.g. ion plasma frequency, ion-sphere radius, coupling parameter, etc.).
                                        
        Parameters
        ----------
        Z : int
            Atomic number 
        ni : float
            Number density [1/cc]
        mi : float
            Atomic mass [g]       
        """

        MI = zbar.MeanIonization(mi, mi*ni, self.T, Z)

        self.Zbar = MI.tf_zbar()
        
        self.ne = self.Zbar * ni # Compute electron number density [1/cm^3]

        self.ai = (4 * np.pi * ni/3)**(-1/3) # Wigner-Seitz radius [cm]
        
        self.gamma = self.Zbar**2 * self.e_squared/(self.ai * self.T) # Coulomb coupling parameter

        self.wp = np.sqrt(4 * np.pi * self.Zbar**2 * ni * self.e_squared/ mi * 1.60218e-12) # Plasma frequency [1/s]

        self.Ef = self.hbar**2/(2 * self.me) * (3 * np.pi**2 * self.ne)**(2/3) * 1.60218e-12 # fermi-energy [ev]

        self.lam_sq = ( self.T**(9/5) + (2/3*self.Ef)**(9/5) )**(5/9) / (4 * np.pi * self.ne * self.e_squared) # [cm]

        self.kappa = self.ai/np.sqrt(self.lam_sq) # inverse screening length # [1]

        self.g  = np.array(self.gamma * np.sqrt( self.kappa**2  + (3*self.gamma)/(1 + 3 * self.gamma) )) # Plasma parameter
        
        
    def knm(self, g, n, m):
        """ Computes the plasma parameters (e.g. ion plasma frequency, ion-sphere radius, coupling parameter, etc.).
                                        
        Parameters
        ----------
        g : float or array_like
            Plasma parameter (eq. 54 from [1]) 
        n : int
            Subscript for collision intergral Knm (eq. C22 from [1])
        m : int
            Subscript for collision integral Knm (eq. C22 from [1])
            
        Returns
        -------
        knm : array_like
            Fit to collision integral (eqs. C22-C24 from [1])
        """

        if n and m == 1:
            a = np.array([1.4660, -1.7836, 1.4313, -0.55833, 0.061162])
            b = np.array([0.081033, -0.091336, 0.051760, -0.50026, 0.17044])
            
        if n and m == 2:
            a = np.array([0.85401, -0.22898, -0.60059, 0.80591, -0.30555])
            b = np.array([0.43475, -0.21147, 0.11116, 0.19665, 0.15195])
            
        knm = np.zeros(len(self.g))
        
        for i in range(len(g)):
            g_arr = np.array([g[i], g[i]**2, g[i]**3, g[i]**4, g[i]**5])
        
            if g[i] < 1:
                knm[i] = -n/4 * np.math.factorial(m - 1) * np.log( np.dot(a,g_arr) ) 
            else:
                knm[i] = (b[0] + b[1]*np.log(g[i]) + b[2]*np.log(g[i])**2)/(1 + b[3]*g[i] + b[4]*g[i]**2)

        return knm

    
    def self_diffusion(self):
        """ Computes the self-diffusion coefficient using eq. 56 from [1].

        Returns
        -------
        D : float or array_like
            Self-diffusion coefficients for system parameters.
            
        Note
        ----
        The structure in which D is returned is in the form of a 3D array (1D if single value inputs). 
        The first axis corresponds to density, the second is temperature, and the third is atomic number. 
        """

        # Single density, single element
        if len(self.num_density.shape) == 1 and len(self.Z) == 1:
            self.plasma_params(self.Z, self.num_density, self.mi)
 
            if self.units_out == 'star':
                D = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1))

            elif self.units_out == 'cgs':
                D = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1)) \
                * self.wp * self.ai**2 
                
            elif self.units_out == 'mks':
                D = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1)) \
                * self.wp * (self.ai/100)**2
                
            else:
                print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')
        
        # Single mass density, multiple elements
        elif len(self.num_density.shape) == 1 and len(self.Z) != 1:
            D = np.zeros([1, len(self.T), len(self.Z)])
            
            for k in range(len(self.Z)):
                self.plasma_params(self.Z[k], self.num_density[k], self.mi[k])

                if self.units_out == 'star':
                    D[0,:,k] = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1))

                elif self.units_out == 'cgs':
                    D[0,:,k] = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1)) \
                    * self.wp * self.ai**2 

                elif self.units_out == 'mks':
                    D[0,:,k] = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1)) \
                    * self.wp * (self.ai/100)**2 

                else:
                    print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')       

        # Multiple mass densities, multiple elements
        else:
            D = np.zeros([self.num_density.shape[1], len(self.T), len(self.Z)])
            
            for k in range(len(self.Z)):
                for i in range(self.num_density.shape[1]):
                    self.plasma_params(self.Z[k], self.num_density[k,i], self.mi[k])

                    if self.units_out == 'star':
                        D[i,:,k] = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1))

                    elif self.units_out == 'cgs':
                        D[i,:,k] = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1))\
                        * self.wp * self.ai**2 
                        
                    elif self.units_out == 'mks':
                        D[i,:,k] = np.sqrt(3*np.pi)/(12 * self.gamma**(5/2) * self.knm(self.g, 1, 1)) \
                        * self.wp * (self.ai/100)**2 
                        
                    else:
                        print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')            
        return D
    
    def viscosity(self):
        """
        Computes the viscosity coefficient using eq. 75 from [1].

        Returns
        -------
        eta : float or array_like
            Viscosity coefficients for system parameters.
            
        Note
        ----
        The structure in which eta is returned is in the form of a 3D array (1D if single value inputs). 
        The first axis corresponds to density, the second is temperature, and the third is atomic number. 
        """
        
        # Single density, single element                              
        if len(self.num_density.shape) == 1 and len(self.Z) == 1:
                           
            self.plasma_params(self.Z, self.num_density, self.mi)
            
            if self.units_out == 'star':
                eta = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) )      

            elif self.units_out == 'cgs':
                eta = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                * self.mi * self.num_density * self.wp * self.ai**2 
                
            elif self.units_out == 'mks':
                eta = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                * self.mi/1000 * (self.num_density/100) * self.wp * (self.ai/100)**2 

            else:
                print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')
                
        # Single mass density, multiple elements        
        elif len(self.num_density.shape) == 1 and len(self.Z) != 1:

            eta = np.zeros([1, len(self.T), len(self.Z)])
            
            for k in range(len(self.Z)):
                self.plasma_params(self.Z[k], self.num_density[k], self.mi[k])

                if self.units_out == 'star':
                    eta[0,:,k] = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) )      

                elif self.units_out == 'cgs':
                    eta[0,:,k] = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                    * self.mi[k] * self.num_density[k] * self.wp * self.ai**2
                   
                elif self.units_out == 'mks':
                    eta[0,:,k] = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                    * self.mi[k]/1000 * (self.num_density[k,i]/100) * self.wp * (self.ai/100)**2
                    
                else:
                    print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"') 
        
        # Multiple mass densities, multiple elements            
        else:
                         
            eta = np.zeros([self.num_density.shape[1], len(self.T), len(self.Z)])

            for k in range(len(self.Z)):
                for i in range(self.num_density.shape[1]):
                    self.plasma_params(self.Z[k], self.num_density[k,i], self.mi[k])

                    if self.units_out == 'star':
                        eta[i,:,k] = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) )      

                    elif self.units_out == 'cgs':
                        eta[i,:,k] = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                        * self.mi[k] * self.num_density[k,i] * self.wp * self.ai**2 
                
                    elif self.units_out == 'mks':
                        eta[i,:,k] = ( 5*np.sqrt(3*np.pi) )/( 36*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                        * self.mi[k]/1000 * (self.num_density[k,i]/100) * self.wp * (self.ai/100)**2
                    else:
                        print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')     
        return eta
        
    def thermal_conductivity(self):
        """
        Computes the ion thermal conductivity coefficient using eq. 82 from [1].

        Returns
        -------
        K : float or array_like
            Ion thermal conductivity coefficients for system parameters.
            
        Note
        ----
        The structure in which K is returned is in the form of a 3D array (1D if single value inputs). 
        The first axis corresponds to density, the second is temperature, and the third is atomic number. 
        """
        
        # Single density, single element
        if len(self.num_density.shape) == 1 and len(self.Z) == 1:
            self.plasma_params(self.Z, self.num_density, self.mi)

            if self.units_out == 'star':
                K = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) )

            elif self.units_out == 'cgs':
                K = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                * self.num_density * self.wp * self.ai**2 * 1.380649e-16 # eV K^{-1}

            elif self.units_out == 'mks':
                K = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                * self.num_density/100 * self.wp * (self.ai/100)**2 * 1.380649e-16 # eV K^{-1}
            else:
                print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')       
                    
        # Single mass density, multiple elements
        elif len(self.num_density.shape) == 1 and len(self.Z) != 1: 
            K = np.zeros([1, len(self.T), len(self.Z)])
            
            for k in range(len(self.Z)):
                self.plasma_params(self.Z[k], self.num_density[k], self.mi[k])
                
                if self.units_out == 'star':
                    K[0,:,k] = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) )

                elif self.units_out == 'cgs':
                    K[0,:,k] = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                    * self.num_density[k] * self.wp * self.ai**2 * 1.380649e-16 # eV K^{-1}
                    
                elif self.units_out == 'mks':
                    K[0,:,k] = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                    * self.num_density[k]/100 * self.wp * (self.ai/100)**2 * 1.380649e-16 # eV K^{-1}
                else:
                    print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')       
        
        # Multiple mass densities, multiple elements                                     
        else:                                 
            K = np.zeros([self.num_density.shape[1], len(self.T), len(self.Z)])
            
            for k in range(len(self.Z)):
                for i in range(self.num_density.shape[1]):
                    self.plasma_params(self.Z[k], self.num_density[k,i], self.mi[k])

                    if self.units_out == 'star':
                        K[i,:,k] = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) )

                    elif self.units_out == 'cgs':
                        K[i,:,k] = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                        * self.num_density[k,i] * self.wp * self.ai**2 * 1.380649e-16 # eV K^{-1}
                        
                    elif self.units_out == 'mks':
                        K[i,:,k] = ( 25*np.sqrt(3*np.pi) )/( 48*self.gamma**(5/2)*self.knm(self.g, 2, 2) ) \
                        * self.num_density[k,i]/100 * self.wp * (self.ai/100)**2 * 1.380649e-16 # eV K^{-1}
                    else:
                        print('Please specify a valid unit system for the returned quantities\nCurrent Options: "star", "cgs", "mks"')          
        return K
    
    
    def plot(self, X, Y, xaxis='temeprature'):
        
        """
        Very rough preliminary plotting method.... in development.
    
        """
        
        if xaxis == 'temperature':
            fig, ax = plt.subplots(nrows=1, ncols=len(self.Z), figsize=(30,8))
            for k in range(len(self.Z)):
                for i in range(self.num_density.shape[1]):
                    ax[k].set_title('Z = ' + str(self.Z[k]), fontsize=22)
                    ax[k].loglog(X, Y[i,:,k], linewidth=2, label='$n_i =$ {:.4e}'.format(self.num_density[k,i]))
                    ax[k].set_xlabel('Temperature', fontsize=20)
                    ax[k].tick_params(axis="x", labelsize=18) 
                    ax[k].tick_params(axis="y", labelsize=18) 
                    ax[k].legend(fontsize=18)

        if xaxis == 'density':
            fig, ax = plt.subplots(nrows=1, ncols=len(self.Z), figsize=(30,8))
            for k in range(len(self.Z)):
                for i in range(len(self.T)):
                    ax[k].set_title(str(self.Z[k]), fontsize=22)
                    ax[k].loglog(X, Y[:,i,k], linewidth=2, label='$T =$ {:.4e}'.format(self.T[i]))
                    ax[k].set_xlabel('Density', fontsize=20)
                    
        return fig, ax


class YVM:
    """Compute viscosity from the Murillo viscosity model.
    
    Parameters
    ----------
    Am : float or arrary_like
        Atomic mass of element (or isotope) in units of grams [g].

    mass_density : float or array_like
        Range of mass densities in units of grams per 
        cubic centimeter [g/cc]. 
    T : float or array_like
        Temperature range in units of electron-volts [eV]

    Z : int or arrray_like
        Atomic number for each element
   
    units_out : str
        Unit system for resulting transport coefficient.
        Default: dimensionless "star" units.

    References
    ----------
    .. [2] `Murillo, M. S. "Viscosity estimates for strongly coupled Yukawa systems." Physical Review E 62.3 (2000): 4115.
        <https://doi.org/10.1103/PhysRevE.62.4115>`_
    """

    def __init__(self, Am, mass_density, T, Z, units_out='star'):
        """ Initialize all parameters to be shared across methods.
        """
        for X in [Am, mass_density, T, Z]:
            if str(type(X)) == "<class 'list'>":
                raise TypeError('type "list" not supported')

        # Check type of input and deal with float cases
        self.mi = Am
        if str(type(self.mi)) != "<class 'numpy.ndarray'>":
            self.mi = np.array([Am])

        self.Z = Z
        if str(type(self.Z)) != "<class 'numpy.ndarray'>":
            self.Z = np.array([Z])
       
        if str(type(mass_density)) != "<class 'numpy.ndarray'>":
            self.num_density = np.array([mass_density])/Am
            print
        else:
            if len(mass_density) > 1:
                self.num_density = np.tile(mass_density, len(Am))/np.repeat(Am, len(mass_density))
                self.num_density = np.reshape(self.num_density, (len(Am) ,len(mass_density)))
            else:
                self.num_density = mass_density/Am
                
        self.T = T
        if str(type(self.T)) != "<class 'numpy.ndarray'>":
            self.T = np.array([T])

        # Check input, throw errors
        if len(self.mi) != len(self.Z):
            raise DimensionMismatchError('atomic mass and nuclear charge arrays must be the same size')

    def viscosity(self):

        # Single density, single element                              
        if len(self.num_density.shape) == 1 and len(self.Z) == 1:

            p = parameters.Params(self.mi, self.mi*self.num_density, self.T, self.Z)
            kappa = p.kappa()
            gamma = p.gamma()
            ifreq = p.wp()
            ai = p.aws()
                           
            A = 0.46*kappa**4/(1 + 0.44*kappa**4)
            B = 1.01*np.exp(-0.92*kappa)
            C = -3.7e-5 + 9e-4*kappa - 2.9e-4*kappa**2

            gamma_ocp = np.abs(A) + np.abs(B)*gamma + np.abs(C)*gamma**2

            lam = 4*np.pi/3 * (3 * gamma_ocp)**(3/2)

            I1 = (180*gamma_ocp*np.pi**(3/2))**(-1)
            I2 = (0.49 - 2.23*gamma_ocp**(-1/3))/(60*np.pi**2)
            I3 = 0.241*gamma_ocp**(1/9)/np.pi**(3/2)

            eta = lam*I1 + (1+lam*I2)**2/(lam*I3)  
            eta *=  self.mi*self.num_density * ifreq * ai**2

                
        # Single mass density, multiple elements        
        elif len(self.num_density.shape) == 1 and len(self.Z) != 1:
            eta = np.zeros([1, len(self.T), len(self.Z)])
            
            for k in range(len(self.Z)):

                p = parameters.Params(self.mi[k], self.mi[k]*self.num_density[k], self.T, self.Z[k])
                kappa = p.kappa()
                gamma = p.gamma()
                ifreq = p.wp()
                ai = p.aws()


                A = 0.46*kappa**4/(1 + 0.44*kappa**4)
                B = 1.01*np.exp(-0.92*kappa)
                C = -3.7e-5 + 9e-4*kappa - 2.9e-4*kappa**2

                gamma_ocp = A + B*gamma + C*gamma**2

                lam = 4*np.pi/3 * (3 * gamma_ocp)**(3/2)

                I1 = (180*gamma_ocp*np.pi**(3/2))**(-1)
                I2 = (0.49 - 2.23*gamma_ocp**(-1/3))/(60*np.pi**2)
                I3 = 0.241*gamma_ocp**(1/9)/np.pi**(3/2)

                eta_val = lam*I1 + (1+lam*I2)**2/(lam*I3)  

                eta[0,:,k] = eta_val * self.mi[k]*self.num_density[k] * ifreq * ai**2
        
        # Multiple mass densities, multiple elements            
        else:
                 
            eta = np.zeros([self.num_density.shape[1], len(self.T), len(self.Z)])
        
            for k in range(len(self.Z)):
                for i in range(self.num_density.shape[1]):

                    p = parameters.Params(self.mi[k], self.mi[k]*self.num_density[k,i], self.T, self.Z[k])
                    kappa = p.kappa()
                    gamma = p.gamma()
                    ifreq = p.wp()
                    ai = p.aws()

                    A = 0.46*kappa**4/(1 + 0.44*kappa**4)
                    B = 1.01*np.exp(-0.92*kappa)
                    C = -3.7e-5 + 9e-4*kappa - 2.9e-4*kappa**2

                    gamma_ocp = A + B*gamma + C*gamma**2

                    lam = 4*np.pi/3 * (3 * gamma_ocp)**(3/2)

                    I1 = (180*gamma_ocp*np.pi**(3/2))**(-1)
                    I2 = (0.49 - 2.23*gamma_ocp**(-1/3))/(60*np.pi**2)
                    I3 = 0.241*gamma_ocp**(1/9)/np.pi**(3/2)

                    eta_val = lam*I1 + (1+lam*I2)**2/(lam*I3)  

                    eta[i,:,k] = eta_val * self.mi[k] * self.num_density[k,i] * ifreq * ai**2
        print('viscosity in units: [g/cm s]')
        return eta


# class MultiSpeciesTrans:

#     def __init__(self, Am, mass_densities, Ts, Zs, units_out='star'):

#         # 0th entry of each following array contains element 1 info
#         # 1st entry of each array contains element 2 info
#         self.mis = Am
#         self.Zs = Zs
#         self.num_densities = mass_densities/Am
#         self.Ts = Ts

#         self.e_squared = 1.4399644800e-7 # [ev cm]
#         self.hbar = 6.5822958e-16 # reduced Planck's constant [eV*s]
#         self.me = 9.1095e-28 # mass of electron [g]
#         self.units_out = units_out

#     def inter_diffusion(self):

#         n = np.sum(self.num_densitites) # Total number density
#         rho = np.sum(self.mass_densities) # Total mass density
#         muij = self.mis[0]*self.mis[1] / np.sum(self.mis) # Reduce mass

#         xi = self.num_densitites/n
#         zi = self.mass_densities/rho

#         atot = (3/(4 * np.pi *n))**(1/3)

#         ne = self.TF_Zbar(self, self.Zs[0], self.num_densitites[0]) * self.num_densitites[0]

#         Ef = self.hbar**2/(2 * self.me) * (3 * np.pi**2 * self.ne)**(2/3) * 1.60218e-12 # fermi-energy [ev]
#         lam_e = ( self.Ts**(9/5) + (2/3*Ef)**(9/5) )**(5/9) / (4 * np.pi * ne * self.e_squared) # [cm]

#         gamm_11 = self.Zs[0]*self.Zs[0]*e_squared / (atot * self.Ts)
#         gamm_22 = self.Zs[2]*self.Zs[2]*e_squared / (atot * self.Ts)
#         gamm_12 = self.Zs[0]*self.Zs[2]*e_squared / (atot * self.Ts)

#         kappa = atot/lam_e

#         g = gam_12 * np.sqrt(kappa + 3*xi[0]**(-1)*gamma_11/(1 + 3*(xi[0]/zi[0])**(1/3) * gamma_11) + \
#                                    + 3*xi[1]**(-1)*gamma_22/(1 + 3*(xi[1]/zi[1])**(1/3) * gamma_22))

#         Dij = 3 * self.Ts**(5/2) / (16 * np.sqrt(2 * np.pi * muij) * n * Zi**2 * Zj**2 * e_squared**2 * Knm(g, 1, 1))
        
#     def tf_zbar(self, Z, ni):

#         """ Thomas Fermi Zbar model.
        
#         Parameters
#         ----------
#         Z : int
#             Atomic number.
            
#         ni : float
#             Number density in units particles per 
#             cubic centimeter [N/cc].
      
#         References
#         ----------
#         Finite Temperature Thomas Fermi Charge State using 
#         R.M. More, "Pressure Ionization, Resonances, and the
#         Continuity of Bound and Free States", Adv. in Atomic 
#         Mol. Phys., Vol. 21, p. 332 (Table IV).
#         """
        
#         alpha = 14.3139
#         beta = 0.6624
#         a1 = 0.003323
#         a2 = 0.9718
#         a3 = 9.26148e-5
#         a4 = 3.10165
#         b0 = -1.7630
#         b1 = 1.43175
#         b2 = 0.31546
#         c1 = -0.366667
#         c2 = 0.983333

#         convert = ni * 1.6726219e-24
#         R = convert/Z
#         T0 = self.T/Z**(4./3.)
#         Tf = T0/(1 + T0)
#         A = a1*T0**a2 + a3*T0**a4
#         B = -np.exp(b0 + b1*Tf + b2*Tf**7)
#         C = c1*Tf + c2
#         Q1 = A*R**B
#         Q = (R**C + Q1**C)**(1/C)
#         x = alpha*Q**beta
        
#         self.Zbar = Z * x/(1 + x + np.sqrt(1 + 2*x))


#     def knm(self, g, n, m):
#         """ Computes the plasma parameters (e.g. ion plasma frequency, ion-sphere radius, coupling parameter, etc.).
                                        
#         Parameters
#         ----------
#         g : float or array_like
#             Plasma parameter (eq. 54 from [1]) 
#         n : int
#             Subscript for collision intergral Knm (eq. C22 from [1])
#         m : int
#             Subscript for collision integral Knm (eq. C22 from [1])
            
#         Returns
#         -------
#         Knm : array_like
#             Fit to collision integral (eqs. C22-C24 from [1])
#         """

#         if n and m == 1:
#             a = np.array([1.4660, -1.7836, 1.4313, -0.55833, 0.061162])
#             b = np.array([0.081033, -0.091336, 0.051760, -0.50026, 0.17044])
            
#         if n and m == 2:
#             a = np.array([0.85401, -0.22898, -0.60059, 0.80591, -0.30555])
#             b = np.array([0.43475, -0.21147, 0.11116, 0.19665, 0.15195])
            
#         Knm = np.zeros(len(self.g))
        
#         for i in range(len(g)):
#             g_arr = np.array([g[i], g[i]**2, g[i]**3, g[i]**4, g[i]**5])
        
#             if g[i] < 1:
#                 Knm[i] = -n/4 * np.math.factorial(m - 1) * np.log( np.dot(a,g_arr) ) 
#             else:
#                 Knm[i] = (b[0] + b[1]*np.log(g[i]) + b[2]*np.log(g[i])**2)/(1 + b[3]*g[i] + b[4]*g[i]**2)

#         return Knm