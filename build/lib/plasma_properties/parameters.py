import numpy as np
import zbar

class Params:


    def __init__(self, Am, mass_density, T, Z, units_out='star'):
        """ Initialize all parameters to be shared across methods.
        """
        
        # Check type of input and deal with float cases
        self.mi = Am
        if str(type(self.mi)) != "<class 'numpy.ndarray'>":
            self.mi = np.array([Am])

        self.Z = Z
        if str(type(self.Z)) != "<class 'numpy.ndarray'>":
            self.Z = np.array([Z])
            
        try:
            self.num_density = np.tile(mass_density, len(Am))/np.repeat(Am, len(mass_density))            
            self.num_density = np.reshape(self.num_density, (len(Am) ,len(mass_density)))
        except:
            self.num_density = np.array([mass_density])/Am

        self.T = T
        if str(type(self.T)) != "<class 'numpy.ndarray'>":
            self.T = np.array([T])

        # Class wide parameter definition in cgs units
        self.e_squared = 1.4399644800e-7 # [ev cm]
        self.hbar = 6.5822958e-16 # reduced Planck's constant [eV*s]
        self.me = 9.1095e-28 # mass of electron [g]
        self.units_out = units_out

    def gamma(self):

        MI = zbar.MeanIonization(self.mi, self.mi*self.num_density, self.T, self.Z)

        Zbar = MI.tf_zbar()
        
        ne = Zbar * self.num_density # Compute electron number density [1/cm^3]

        ai = (4 * np.pi * self.num_density/3)**(-1/3) # Wigner-Seitz radius [cm]
        
        g = Zbar**2 * self.e_squared/(ai * self.T) # Coulomb coupling parameter

        return g