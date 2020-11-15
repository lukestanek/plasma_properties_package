import numpy as np
import zbar

class Params:
    """Compute common plasma parameters.
    
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
    """


    def __init__(self, Am, mass_density, T, Z):
        """ Initialize all parameters to be shared across methods.
        """
        
        # Check type of input and deal with float cases
        self.mi = Am
        self.mass_density = mass_density
        self.num_density  = self.mass_density/self.mi
        self.T = T
        self.Z = Z    

        # Class wide parameter definition in cgs units
        self.e_squared = 1.4399644800e-7 # [ev cm]
        self.hbar = 6.5822958e-16 # reduced Planck's constant [eV*s]
        self.me = 9.1095e-28 # mass of electron [g]
        # self.units_out = units_out

    def wp(self):
        """ Ion plasma frequency

        .. math::

            \omega_p = \\left(\\dfrac{4  n_i  \\langle Z \\rangle^2 e^2}{m_i}\\right)^{1/2} 
        """


        MI = zbar.MeanIonization(self.mi, self.mass_density, self.T, self.Z)

        Zbar = MI.tf_zbar()

        # Compute ion plasma frequency
        wpi = ((4*np.pi * self.num_density * Zbar**2 * self.e_squared)/self.mi * 1.602e-12 )**0.5 # 1/s

        return wpi

    def aws(self):
        """ Wigner-Seitz radius

        .. math::

            a_{ws}= (4  \pi n_i/3)^{-1/3} 
        """
        return (4 * np.pi * self.num_density/3)**(-1/3) # cm

    def degen(self):
        """ degeneracy  parameter

        .. math::

            \\theta = E_f/T

        where :math:`E_F` is  the Fermi energy :math:`E_F = \\frac{\hbar^2}{2m_e}(3 \pi^2 n_e)^{2/3}`.
        """

        MI = zbar.MeanIonization(self.mi, self.mass_density, self.T, self.Z)

        Zbar = MI.tf_zbar()

        ne = Zbar*self.num_density

        Ef = self.hbar**2/(2*self.me)*(3*np.pi**2*ne)**(2/3) * 1.60218e-12
        theta = Ef/self.T

        return theta

    def gamma(self):
        """ Coulomb coupling parameter

        .. math::

            \\Gamma = \\dfrac{\\langle Z \\rangle^2 e^2}{a_i T}
        """
        MI = zbar.MeanIonization(self.mi, self.mass_density, self.T, self.Z)

        Zbar = MI.tf_zbar()
        
        ne = Zbar * self.num_density # Compute electron number density [1/cm^3]

        ai = (4 * np.pi * self.num_density/3)**(-1/3) # Wigner-Seitz radius [cm]
        
        g = Zbar**2 * self.e_squared/(ai * self.T) # Coulomb coupling parameter

        return g

    def kappa(self):
        """ Inverse electron screening length

        .. math::

            \\kappa = \\dfrac{a_i}{\\lambda_{TF}},

        where :math:`\\lambda_{TF}` is the Thomas-Fermi screening length.
        """

        MI = zbar.MeanIonization(self.mi, self.mass_density, self.T, self.Z)

        Zbar = MI.tf_zbar()
        
        ne = Zbar * self.num_density # Compute electron number density [1/cm^3]

        ai = (4 * np.pi * self.num_density/3)**(-1/3) # Wigner-Seitz radius [cm]

        Ef = self.hbar**2/(2 * self.me) * (3 * np.pi**2 * ne)**(2/3) * 1.60218e-12 # fermi-energy [ev]

        lam_sq = ( self.T**(9/5) + (2/3*Ef)**(9/5) )**(5/9) / (4 * np.pi * ne * self.e_squared) # [cm]

        k = ai/np.sqrt(lam_sq) # inverse screening length # [1]

        return k

