import numpy as np

# def ZBAR(Am, mass_density, T, Z, model):
#     """ Compute the mean ionization state (<Z> or Zbar) for a given system.

#     Parameters
#     ----------
#     Am : float or arrary_like
#         Atomic mass of element (or isotope) in units of grams [g].
        
#     mass_density : float or array_like
#         Range of mass densities in units of grams per 
#         cubic centimeter [g/cc]. 
        
#     T : float or array_like
#         Temperature range in units of elevtron-volts [eV]
        
#     Z : int or arrray_like
#         Atomic number for each element
#     """

#     # Instantiate the MeanIonization object
#     MI = MeanIonization(Am, mass_density, T, Z, model)

#     # Catch incorrect models
#     try:
#         return MI.get_zbar()
#     except:
#         print('Error in Choice of Model')

class MeanIonization:
    """ Compute the mean ionization state (<Z> or Zbar) for a given system.

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
    """
    def __init__(self, Am, mass_density, T, Z):
        
        # Check type of input and deal with float cases
        if str(type(Am)) != "<class 'numpy.ndarray'>":
            self.Am = np.array([Am])
        else:
            self.Am = Am

        if str(type(Z)) != "<class 'numpy.ndarray'>":
            self.Z = np.array([Z])
        else:
            self.Z = Z
            
        if str(type(mass_density)) != "<class 'numpy.ndarray'>":
            self.num_density = np.array([mass_density])/Am
        else:
            if len(mass_density) > 1:
                self.num_density = np.tile(mass_density, len(Am))/np.repeat(Am, len(mass_density))
                self.num_density = np.reshape(self.num_density, (len(Am) ,len(mass_density)))
            else:
                self.num_density = mass_density/Am

        if str(type(T)) != "<class 'numpy.ndarray'>":
            self.T = np.array([T])
        else:
            self.T = T

        # # Select model for computing Zbar
        # if model=='TF':
        #     self.zbar = self.tf_zbar()
        # elif model=='SAHA':
        #     print("Coming Soon")

    # def get_zbar(self):
    #     return self.zbar

    def tf_zbar(self):
            """ Thomas Fermi Zbar model.

            References
            ----------
            Finite Temperature Thomas Fermi Charge State using 
            R.M. More, "Pressure Ionization, Resonances, and the
            Continuity of Bound and Free States", Adv. in Atomic 
            Mol. Phys., Vol. 21, p. 332 (Table IV).
            """
            
            # Single mass density, single element
            if len(self.num_density.shape) == 1 and len(self.Z) == 1:
                
                alpha = 14.3139
                beta = 0.6624
                a1 = 0.003323
                a2 = 0.9718
                a3 = 9.26148e-5
                a4 = 3.10165
                b0 = -1.7630
                b1 = 1.43175
                b2 = 0.31546
                c1 = -0.366667
                c2 = 0.983333

                convert = self.num_density * 1.6726219e-24
                R = convert/self.Z
                T0 = self.T/self.Z**(4./3.)
                Tf = T0/(1 + T0)
                A = a1*T0**a2 + a3*T0**a4
                B = -np.exp(b0 + b1*Tf + b2*Tf**7)
                C = c1*Tf + c2
                Q1 = A*R**B
                Q = (R**C + Q1**C)**(1/C)
                x = alpha*Q**beta

                Zbar = self.Z * x/(1 + x + np.sqrt(1 + 2*x))
        
            # Single mass density, multiple elements
            elif len(self.num_density.shape) == 1 and len(self.Z) != 1:

                Zbar = np.zeros([1, len(self.T), len(self.Z)])
               
                for k in range(len(self.Z)):
                    alpha = 14.3139
                    beta = 0.6624
                    a1 = 0.003323
                    a2 = 0.9718
                    a3 = 9.26148e-5
                    a4 = 3.10165
                    b0 = -1.7630
                    b1 = 1.43175
                    b2 = 0.31546
                    c1 = -0.366667
                    c2 = 0.983333

                    convert = self.num_density[k] * 1.6726219e-24
                    R = convert/self.Z[k]
                    T0 = self.T/self.Z[k]**(4./3.)
                    Tf = T0/(1 + T0)
                    A = a1*T0**a2 + a3*T0**a4
                    B = -np.exp(b0 + b1*Tf + b2*Tf**7)
                    C = c1*Tf + c2
                    Q1 = A*R**B
                    Q = (R**C + Q1**C)**(1/C)
                    x = alpha*Q**beta

                    Zbar[0,:,k] = self.Z[k] * x/(1 + x + np.sqrt(1 + 2*x))

            # Multiple mass densities, multiple elements
            else:
                Zbar = np.zeros([self.num_density.shape[1], len(self.T), len(self.Z)])
                
                for k in range(len(self.Z)):
                    for i in range(self.num_density.shape[1]):
                        
                        alpha = 14.3139
                        beta = 0.6624
                        a1 = 0.003323
                        a2 = 0.9718
                        a3 = 9.26148e-5
                        a4 = 3.10165
                        b0 = -1.7630
                        b1 = 1.43175
                        b2 = 0.31546
                        c1 = -0.366667
                        c2 = 0.983333

                        convert = self.num_density[k,i] * 1.6726219e-24
                        R = convert/self.Z[k]
                        T0 = self.T/self.Z[k]**(4./3.)
                        Tf = T0/(1 + T0)
                        A = a1*T0**a2 + a3*T0**a4
                        B = -np.exp(b0 + b1*Tf + b2*Tf**7)
                        C = c1*Tf + c2
                        Q1 = A*R**B
                        Q = (R**C + Q1**C)**(1/C)
                        x = alpha*Q**beta

                        Zbar[i,:,k] = self.Z[k] * x/(1 + x + np.sqrt(1 + 2*x))

            return Zbar
