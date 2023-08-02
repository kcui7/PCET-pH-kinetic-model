import numpy as np

class GeneralizedSquareDiagram(object):
    """
    This class stores all thermodynamic parameters which will be used to calculate the mole fractions at given pH. 
    This class also enables the use of potential-dependent parameters by
    using callable functions as the K or E_ET1 parameters. 

    Input description and data type:
    KaOx (float or function): acid dissociation equilibrium constant of the oxidized species 
    KaRed (float or function): acid dissociation equilibrium constant of the reduced species 
    K1 (float or function): association equilibrium constant of the oxidized species with H3O+ 
    K2 (float or function): association equilibrium constant of the reduced species with OH-
    E_ET1 (float or function): redox potential of the ET1 channel (in Volts)
    T (float): the temperature, default value = 298K
    Kw (float): water self-dissociation equilibrium constant at the given temperature, default value = 1e-14
    """
    def __init__(self, KaOx=1e-3, KaRed=1e-9, K1=1e2, K2=1e3, E_ET1=1.0, T=298, Kw=1e-14):
        self.T = T
        self.R = 8.314
        self.F = 96485
        self.parameters = {'KaOx':KaOx,
                           'KaRed':KaRed,
                           'K1':K1,
                           'K2':K2,
                           'Kw':Kw,
                           'E_ET1':E_ET1,
                           }

        self.pot_dep_parameters = {key:self.parameters[key] for key in self.parameters.keys() if callable(self.parameters[key])} 
        self.pot_indep_parameters = {key:self.parameters[key] for key in self.parameters.keys() if not callable(self.parameters[key])} 

        # When use potential-dependent parameters
        # During initialization, KaOx and E_ET1 must be both potential-dependent to ensure 
        # the free energy of the overall reaction: Ox + H+ + e- <==> RedH is potential-independent. 
        # On the other hand, for the other legs, say if KaRed is set to be potential-dependent, 
        # the code will automatically make E_ET2 also potential-dependent. 

        if not (callable(self.parameters['E_ET1']) == callable(self.parameters['KaOx'])): 
            raise RuntimeError("KaOx and E_ET1 must be both potential-dependent to ensure the standard free energy of the overall reaction: Ox + H+ + e- <==> RedH is potential-independent. ")


        # detect if any functions are used as parameters, if so, enable potential-dependent parameters
        if len(self.pot_dep_parameters.keys()) == 0:
            self.potential_dependence = False
        else:
            self.potential_dependence = True


        if self.potential_dependence == False:
            E_ET2 = self.parameters['E_ET1'] - self.R*self.T/self.F * np.log(self.parameters['KaOx'] / self.parameters['KaRed'])
            E_EPT1 = self.parameters['E_ET1'] - self.R*self.T/self.F * np.log(self.parameters['KaOx'] * self.parameters['K1'])
            E_EPT2 = E_ET2 - self.R*self.T/self.F * np.log(self.parameters['KaRed'] / self.parameters['K2'] / self.parameters['Kw'])
            self.parameters['E_ET2'] = E_ET2 
            self.parameters['E_EPT1'] = E_EPT1 
            self.parameters['E_EPT2'] = E_EPT2
            
        else:
            # transfer potential-independent parameters to constant functions
            if not callable(self.parameters['KaOx']):
                def KaOx_func(E):
                    return self.pot_indep_parameters['KaOx']
                self.parameters['KaOx'] = KaOx_func
            if not callable(self.parameters['KaRed']):
                def KaRed_func(E):
                    return self.pot_indep_parameters['KaRed']
                self.parameters['KaRed'] = KaRed_func
            if not callable(self.parameters['K1']):
                def K1_func(E):
                    return self.pot_indep_parameters['K1']
                self.parameters['K1'] = K1_func
            if not callable(self.parameters['K2']):
                def K2_func(E):
                    return self.pot_indep_parameters['K2']
                self.parameters['K2'] = K2_func
            if not callable(self.parameters['Kw']):
                def Kw_func(E):
                    return self.pot_indep_parameters['Kw']
                self.parameters['Kw'] = Kw_func
            if not callable(self.parameters['E_ET1']):
                def E_ET1_func(E):
                    return self.pot_indep_parameters['E_ET1']
                self.parameters['E_ET1'] = E_ET1_func

            def E_ET2_func(E):
                return self.parameters['E_ET1'](E) - self.R*self.T/self.F * np.log(self.parameters['KaOx'](E) / self.parameters['KaRed'](E))
            def E_EPT1_func(E):
                return self.parameters['E_ET1'](E) - self.R*self.T/self.F * np.log(self.parameters['KaOx'](E) * self.parameters['K1'](E))
            def E_EPT2_func(E):
                return E_ET2_func(E) - self.R*self.T/self.F * np.log(self.parameters['KaRed'](E) / self.parameters['K2'](E) / self.parameters['Kw'](E))

            self.parameters['E_ET2'] = E_ET2_func 
            self.parameters['E_EPT1'] = E_EPT1_func 
            self.parameters['E_EPT2'] = E_EPT2_func


    def __repr__(self):
        if self.potential_dependence == False:
            str = "The following parameters are defined in the generlized square diagram\n"
            for parameter in ['KaOx', 'KaRed', 'K1', 'K2', 'Kw', 'E_ET1', 'E_ET2', 'E_EPT1', 'E_EPT2']:
                if parameter.startswith('K'):
                    str += "%s = %.2e\n"%(parameter, self.parameters[parameter])
                elif parameter.startswith('E'):
                    str += "%s = %.3f V\n"%(parameter, self.parameters[parameter])
                else:
                    pass

            return str
        else:
            str = "The following functions are defined for the potential-dependent parameters in the generlized square diagram\n"
            for parameter in ['KaOx', 'KaRed', 'K1', 'K2', 'Kw', 'E_ET1', 'E_ET2', 'E_EPT1', 'E_EPT2']:
                str += "%s = %s\n"%(parameter, self.parameters[parameter].__repr__())

            return str


    def get_parameters(self, key='all'):
        if key == 'all':
            return self.parameters
        else:
            return self.parameters[key]


    def set_parameters(self, **kwargs): 
        # re-initialize the instance
        # first, remove E_ET2, E_EPT1, and E_EPT2 from self.parameters
        # they will be re calculated later

        del self.parameters['E_ET2']    
        del self.parameters['E_EPT1']    
        del self.parameters['E_EPT2']    

        parameters = {key:self.parameters[key] for key in self.parameters.keys()}
        T = self.T

        # second, update parameters from **kwargs 
        for key in kwargs.keys():
            if key in self.parameters.keys() and key.startswith('K'):
                parameters[key] = kwargs[key]
            elif key == 'E_ET1':
                parameters[key] = kwargs[key]
            elif key in ('E_ET2', 'E_EPT1', 'E_EPT2'):
                print("!!!\nWARNING: '%s' is calculated based on other parameters.\n!!!"%key)
            elif key == 'T':
                T = kwargs[key]
            else:
                pass

        # finally, reinitialize by calling __init__()
        self.__init__(KaOx=parameters['KaOx'],
                      KaRed=parameters['KaRed'],
                      K1=parameters['K1'],
                      K2=parameters['K2'],
                      E_ET1=parameters['E_ET1'],
                      T=T,
                      Kw=parameters['Kw'],
                      )


    def calc_mole_fraction(self, pH, E=None, species='all'):
        """
        This method calculate the mole fraction of the specified species as a function of pH
        With potential-dependent parameters, the user should also specify E (float). 

        Input description:
        pH (float or np.ndarray): the pH of interest
        E (float or np.ndarray): the electrode potential of interest, default value = None
        !!! NOTE: pH and E cannot be both np.ndarray
        species (string): the species of interest, possible choices: 'OxH', 'Ox', 'OxH3O', 'RedH', 'Red', 'RedHOH', 'all'

        returns:
        if species == 'all': mole fractions of all species, with the order of: 'OxH', 'Ox', 'OxH3O', 'RedH', 'Red', 'RedHOH'
        else: mole fraction of the specific species. 
        """
        if isinstance(pH,np.ndarray) and isinstance(E,np.ndarray) and len(pH) != 1 and len(E) !=1:
            raise TypeError("pH and E cannot be both arrays")

        cH = 10**(-pH)

        if self.potential_dependence == False:
            KaOx = self.parameters['KaOx']
            KaRed = self.parameters['KaRed']
            K1 = self.parameters['K1']
            K2 = self.parameters['K2']
            Kw = self.parameters['Kw']
        elif self.potential_dependence == True and not(isinstance(E,float) or isinstance(E,np.ndarray) or isinstance(E,int)):
            raise ValueError("Please specify tht electrode potential. ")
        else: 
            KaOx = self.parameters['KaOx'](E)
            KaRed = self.parameters['KaRed'](E)
            K1 = self.parameters['K1'](E)
            K2 = self.parameters['K2'](E)
            Kw = self.parameters['Kw'](E)

        chi_OxH = 1/(1 + KaOx/cH + K1*KaOx) 
        chi_Ox = 1/(1 + cH/KaOx + K1*cH)
        chi_OxH3O = 1/(1 + 1/KaOx/K1 + 1/K1/cH)
        chi_RedH = 1/(1 + KaRed/cH + K2*Kw/cH)
        chi_Red = 1/(1 + cH/KaRed + K2*Kw/KaRed)
        chi_RedHOH = 1/(1 + cH/K2/Kw + KaRed/K2/Kw)

        if species == 'all':
            return chi_OxH, chi_Ox, chi_OxH3O, chi_RedH, chi_Red, chi_RedHOH
        elif species == 'OxH':
            return chi_OxH
        elif species == 'Ox':
            return chi_Ox
        elif species == 'OxH3O':
            return chi_OxH3O
        elif species == 'RedH':
            return chi_RedH
        elif species == 'Red':
            return chi_Red
        elif species == 'RedHOH':
            return chi_RedHOH
        else:
            raise KeyError("'%s'"%species)


    def calc_apparent_potential(self, pH):
        if self.potential_dependence == False:
            chi_OxH = self.calc_mole_fraction(pH, species='OxH')
            chi_RedH = self.calc_mole_fraction(pH, species='RedH')
            Eapp = self.parameters['E_ET1'] - self.R*self.T/self.F * np.log(chi_RedH/chi_OxH) 

            return Eapp
        else:
            cH = 10**(-pH)
            Eold = 0

            KaOx = self.parameters['KaOx'](Eold)
            KaRed = self.parameters['KaRed'](Eold)
            K1 = self.parameters['K1'](Eold)
            K2 = self.parameters['K2'](Eold)
            Kw = self.parameters['Kw'](Eold)
            chi_OxH = 1/(1 + KaOx/cH + K1*KaOx) 
            chi_RedH = 1/(1 + KaRed/cH + K2*Kw/cH)

            Enew = self.parameters['E_ET1'](Eold) - self.R*self.T/self.F * np.log(chi_RedH/chi_OxH) 

            while np.any(np.abs(Enew-Eold) > 0.001):
                Eold = Enew

                KaOx = self.parameters['KaOx'](Eold)
                KaRed = self.parameters['KaRed'](Eold)
                K1 = self.parameters['K1'](Eold)
                K2 = self.parameters['K2'](Eold)
                Kw = self.parameters['Kw'](Eold)
                chi_OxH = 1/(1 + KaOx/cH + K1*KaOx) 
                chi_RedH = 1/(1 + KaRed/cH + K2*Kw/cH)

                Enew = self.parameters['E_ET1'](Eold) - self.R*self.T/self.F * np.log(chi_RedH/chi_OxH) 

            return Enew

