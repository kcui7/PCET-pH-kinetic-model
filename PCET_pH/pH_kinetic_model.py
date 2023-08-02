import numpy as np
from scipy.special import erf
from .square_diagram import GeneralizedSquareDiagram


class pHKineticModel(GeneralizedSquareDiagram):
    """
    This class stores all thermodynamics and kinetics parameters in the PCET kinetic model which can be used to 
    calculate the mole fractions, apparent potentials and apparent rate constants as functions of pH (and apparent potential). 
    This class can be used to describe both homogeneous and electrochemical reactions. 
    For the later, we pre-defined that the potential-dependent rate constants are described by either Butler-Volmer or Marcus-Gerischer model.
    This class also enables the use of potential-dependent and independent thermodynamics paramters for electrochemical systems

    Input description and data type:
    KaOx (float or function): acid dissociation equilibrium constant of the oxidized species 
    KaRed (float or function): acid dissociation equilibrium constant of the reduced species 
    K1 (float or function): association equilibrium constant of the oxidized species with H3O+ 
    K2 (float or function): association equilibrium constant of the reduced species with OH-
    E_ET1 (float or function): redox potential of the ET1 channel (in Volts)
    T (float): the temperature, default value = 298K
    Kw (float): water self-dissociation equilibrium constant at the given temperature, default value = 1e-14
    
    k_ET1/k_ET2/k_EPT1/k_EPT2 (float or tuple): rate constants (and other parameters) of the ET1/ET2/EPT1/EPT2 channel
                                                if reaction = 'homogeneous': k = rate constant (float)
                                                if reaction = 'BV': k = (rate constant, cathodic charge transfer coefficient)
                                                if reaction = 'MG': k = (rate constant, reorganization energy (in eV))
    kinetics (string): choose among 'homogeneous', 'Butler-Volmer' or 'BV', and 'Marcus-Gerischer' or 'MG'
    reaction (string): choose between 'reduction' and 'oxidation'
    """
    def __init__(self, KaOx=1e-3, KaRed=1e-9, K1=1e2, K2=1e3, E_ET1=1.0, T=298, Kw=1e-14,
                       k_ET1=1e5, k_ET2=1e5, k_EPT1=1e5, k_EPT2=1e5,
                       reaction='reduction', kinetics='homogeneous'):
        super().__init__(KaOx=KaOx, KaRed=KaRed, K1=K1, K2=K2, E_ET1=E_ET1, T=T, Kw=Kw)

        self.rate_constants = {}
        self.reaction = reaction
        self.kinetics = kinetics

        if self.kinetics == 'homogeneous':
            if self.potential_dependence == True:
                raise ValueError("The thermodynamics parameters are potential-dependent. Please choose the reaction kinetics among 'Marcus-Gerischer' ('MG') and 'Butler-Volmer' ('BV'). ")
            if not (isnumber(k_ET1) and isnumber(k_ET2) and isnumber(k_EPT1) and isnumber(k_EPT2)):
                raise TypeError("The system is homogeneous, please input k = rate constant (float). ")
            else:
                self.rate_constants['k_ET1'] = k_ET1
                self.rate_constants['k_ET2'] = k_ET2
                self.rate_constants['k_EPT1'] = k_EPT1
                self.rate_constants['k_EPT2'] = k_EPT2

        elif self.kinetics == 'Marcus-Gerischer' or self.kinetics == 'MG':
            # check parameter format
            if not (isarray(k_ET1) and isarray(k_ET2) and isarray(k_EPT1) and isarray(k_EPT2)):
                raise TypeError("The system is electrochemical with Marcus-Gerischer kinetics, please input k = (rate constant, reorganization energy (in eV) ). ")
            elif not (len(k_ET1)==2 and len(k_ET2)==2 and len(k_EPT1)==2 and len(k_EPT2)==2):
                raise TypeError("The system is electrochemical with Marcus-Gerischer kinetics, please input k = (rate constant, reorganization energy (in eV) ). ")
            else:
                self.max_rate_constants = {}
                self.reorganization_energies = {}

                self.max_rate_constants['k_ET1'], self.reorganization_energies['k_ET1'] = k_ET1
                self.max_rate_constants['k_ET2'], self.reorganization_energies['k_ET2'] = k_ET2
                self.max_rate_constants['k_EPT1'], self.reorganization_energies['k_EPT1'] = k_EPT1
                self.max_rate_constants['k_EPT2'], self.reorganization_energies['k_EPT2'] = k_EPT2

                # define potential-dependent rate constants based on whether the reaction is reduction or oxidation,
                # and the thermodynamic parameters are potential-dependent or not. 
                if self.reaction == 'reduction':
                    if self.potential_dependence == False:
                        def k_ET1_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET1']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET1'] / self.R/self.T

                            k = self.max_rate_constants['k_ET1']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_ET2_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET2']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET2'] / self.R/self.T

                            k = self.max_rate_constants['k_ET2']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT1_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT1'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT1']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT2_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT2'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT2']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k
                    else:
                        def k_ET1_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET1'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET1'] / self.R/self.T

                            k = self.max_rate_constants['k_ET1']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_ET2_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET2'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET2'] / self.R/self.T

                            k = self.max_rate_constants['k_ET2']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT1_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT1'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT1']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT2_MG_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT2'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT2']*0.5*(1.-erf((eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                    self.rate_constants['k_ET1'] = k_ET1_MG_red_func
                    self.rate_constants['k_ET2'] = k_ET2_MG_red_func
                    self.rate_constants['k_EPT1'] = k_EPT1_MG_red_func
                    self.rate_constants['k_EPT2'] = k_EPT2_MG_red_func

                elif self.reaction == 'oxidation':
                    if self.potential_dependence == False:
                        def k_ET1_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET1']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET1'] / self.R/self.T

                            k = self.max_rate_constants['k_ET1']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_ET2_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET2']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET2'] / self.R/self.T

                            k = self.max_rate_constants['k_ET2']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT1_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT1'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT1']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT2_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2']) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT2'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT2']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k
                    else:
                        def k_ET1_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET1'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET1'] / self.R/self.T

                            k = self.max_rate_constants['k_ET1']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_ET2_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET2'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_ET2'] / self.R/self.T

                            k = self.max_rate_constants['k_ET2']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT1_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT1'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT1']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                        def k_EPT2_MG_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2'](E)) / self.R/self.T
                            Lambda = self.F * self.reorganization_energies['k_EPT2'] / self.R/self.T

                            k = self.max_rate_constants['k_EPT2']*0.5*(1.-erf((-eta+Lambda)/np.sqrt(4*Lambda) ))
                            return k

                    self.rate_constants['k_ET1'] = k_ET1_MG_ox_func
                    self.rate_constants['k_ET2'] = k_ET2_MG_ox_func
                    self.rate_constants['k_EPT1'] = k_EPT1_MG_ox_func
                    self.rate_constants['k_EPT2'] = k_EPT2_MG_ox_func

                else:
                    raise ValueError("reaction should be either 'reduction' or 'oxidation'. ")

        elif self.kinetics == 'Butler-Volmer' or self.kinetics == 'BV':
            # check parameter format
            if not (isarray(k_ET1) and isarray(k_ET2) and isarray(k_EPT1) and isarray(k_EPT2)):
                raise TypeError("The system is electrochemical with Butler-Volmer kinetics, please input k = (rate constant, cathodic charge transfer coefficient). ")
            elif not (len(k_ET1)==2 and len(k_ET2)==2 and len(k_EPT1)==2 and len(k_EPT2)==2):
                raise TypeError("The system is electrochemical with Butler-Volmer kinetics, please input k = (rate constant, cathodic charge transfer coefficient). ")
            else:
                self.standard_rate_constants = {}
                self.charge_transfer_coefficients = {}

                self.standard_rate_constants['k_ET1'], self.charge_transfer_coefficients['k_ET1'] = k_ET1
                self.standard_rate_constants['k_ET2'], self.charge_transfer_coefficients['k_ET2'] = k_ET2
                self.standard_rate_constants['k_EPT1'], self.charge_transfer_coefficients['k_EPT1'] = k_EPT1
                self.standard_rate_constants['k_EPT2'], self.charge_transfer_coefficients['k_EPT2'] = k_EPT2

                # define potential-dependent rate constants based on whether the reaction is reduction or oxidation,
                # and the thermodynamic parameters are potential-dependent or not. 
                if self.reaction == 'reduction':
                    if self.potential_dependence == False:
                        def k_ET1_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET1']) / self.R/self.T
                            k = self.standard_rate_constants['k_ET1']*np.exp(-self.charge_transfer_coefficients['k_ET1']*eta)
                            return k

                        def k_ET2_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET2']) / self.R/self.T
                            k = self.standard_rate_constants['k_ET2']*np.exp(-self.charge_transfer_coefficients['k_ET2']*eta)
                            return k

                        def k_EPT1_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1']) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT1']*np.exp(-self.charge_transfer_coefficients['k_EPT1']*eta)
                            return k

                        def k_EPT2_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2']) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT2']*np.exp(-self.charge_transfer_coefficients['k_EPT2']*eta)
                            return k
                    else:
                        def k_ET1_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET1'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_ET1']*np.exp(-self.charge_transfer_coefficients['k_ET1']*eta)
                            return k

                        def k_ET2_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_ET2'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_ET2']*np.exp(-self.charge_transfer_coefficients['k_ET2']*eta)
                            return k

                        def k_EPT1_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT1']*np.exp(-self.charge_transfer_coefficients['k_EPT1']*eta)
                            return k

                        def k_EPT2_BV_red_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT2']*np.exp(-self.charge_transfer_coefficients['k_EPT2']*eta)
                            return k

                    self.rate_constants['k_ET1'] = k_ET1_BV_red_func
                    self.rate_constants['k_ET2'] = k_ET2_BV_red_func
                    self.rate_constants['k_EPT1'] = k_EPT1_BV_red_func
                    self.rate_constants['k_EPT2'] = k_EPT2_BV_red_func

                elif self.reaction == 'oxidation':
                    if self.potential_dependence == False:
                        def k_ET1_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET1']) / self.R/self.T
                            k = self.standard_rate_constants['k_ET1']*np.exp( (1-self.charge_transfer_coefficients['k_ET1']) * eta)
                            return k

                        def k_ET2_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET2']) / self.R/self.T
                            k = self.standard_rate_constants['k_ET2']*np.exp( (1-self.charge_transfer_coefficients['k_ET2']) * eta)
                            return k

                        def k_EPT1_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1']) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT1']*np.exp( (1-self.charge_transfer_coefficients['k_EPT1']) * eta)
                            return k

                        def k_EPT2_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2']) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT2']*np.exp( (1-self.charge_transfer_coefficients['k_EPT2']) * eta)
                            return k
                    else:
                        def k_ET1_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET1'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_ET1']*np.exp( (1-self.charge_transfer_coefficients['k_ET1']) * eta)
                            return k

                        def k_ET2_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_ET2'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_ET2']*np.exp( (1-self.charge_transfer_coefficients['k_ET2']) * eta)
                            return k

                        def k_EPT1_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT1'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT1']*np.exp( (1-self.charge_transfer_coefficients['k_EPT1']) * eta)
                            return k

                        def k_EPT2_BV_ox_func(E):
                            eta = self.F * (E - self.parameters['E_EPT2'](E)) / self.R/self.T
                            k = self.standard_rate_constants['k_EPT2']*np.exp( (1-self.charge_transfer_coefficients['k_EPT2']) * eta)
                            return k

                    self.rate_constants['k_ET1'] = k_ET1_BV_ox_func
                    self.rate_constants['k_ET2'] = k_ET2_BV_ox_func
                    self.rate_constants['k_EPT1'] = k_EPT1_BV_ox_func
                    self.rate_constants['k_EPT2'] = k_EPT2_BV_ox_func

                else:
                    raise ValueError("reaction should be either 'reduction' or 'oxidation'. ")
        else:
            raise NotImplementedError("Please choose the reaction kinetics among 'homogeneous', 'Marcus-Gerischer' ('MG'), and 'Butler-Volmer' ('BV'). ")


    def __repr__(self):
        str = super().__repr__()

        if self.kinetics == 'homogeneous':
            str2 = "\nThe following rate constants are defined for the homogeneous system. \n"
            for key in ['ET1', 'ET2', 'EPT1', 'EPT2']:
                str2 += "%s: %f\n"%(key, self.rate_constants['k_'+key])
        elif self.kinetics == 'Butler-Volmer' or self.kinetics == 'BV':
            str2 = "\nThe following standard rate constants and cathodic charge transfer coefficients are defined for the electrochemical system. \n"
            for key in ['ET1', 'ET2', 'EPT1', 'EPT2']:
                str2 += "%s: %f, %.3f\n"%(key, self.standard_rate_constants['k_'+key], self.charge_transfer_coefficients['k_'+key])
        elif self.kinetics == 'Marcus-Gerischer' or self.kinetics == 'MG':
            str2 = "\nThe following maximum rate constants and reorgnization energies are defined for the electrochemical system. \n"
            for key in ['ET1', 'ET2', 'EPT1', 'EPT2']:
                str2 += "%s: %f, %.3f eV\n"%(key, self.max_rate_constants['k_'+key], self.reorganization_energies['k_'+key])
        else:
            pass

        str += str2

        return str


    def get_parameters(self, key='all'):
        if key == 'all':
            tmp_dict = {}
            for tmp_key in self.parameters.keys():
                tmp_dict[tmp_key] = self.parameters[tmp_key]
            for tmp_key in self.rate_constants.keys():
                tmp_dict[tmp_key] = self.rate_constants[tmp_key]
            return tmp_dict
        if key.startswith('K') or key.startswith('E'):
            return self.parameters[key]
        elif key.startswith('k'):
            return self.rate_constants[key]
        else:
            raise KeyError("'%s'"%key)


    def set_parameters(self, **kwargs):
        # re-initialize the instance
        # first, remove E_ET2, E_EPT1, and E_EPT2 from self.parameters
        # they will be re calculated later

        del self.parameters['E_ET2']    
        del self.parameters['E_EPT1']    
        del self.parameters['E_EPT2']    

        parameters = {key:self.parameters[key] for key in self.parameters.keys()}
        T = self.T
        reaction = self.reaction
        kinetics = self.kinetics
        if kinetics == 'homogeneous':
            rate_constants = {key:self.rate_constants[key] for key in self.rate_constants.keys()}
        elif kinetics == 'Butler-Volmer' or kinetics == 'BV':
            rate_constants = {key:(self.standard_rate_constants[key], self.charge_transfer_coefficients[key]) for key in self.rate_constants.keys()}
        elif kinetics == 'Marcus-Gerischer' or kinetics == 'MG':
            rate_constants = {key:(self.max_rate_constants[key], self.reorganization_energies[key]) for key in self.rate_constants.keys()}
        else:
            pass

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
            elif key in self.rate_constants.keys():
                rate_constants[key] = kwargs[key]
            elif key == 'reaction':
                reaction = kwargs[key]
            elif key == 'kinetics':
                kinetics = kwargs[key]
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
                      k_ET1=rate_constants['k_ET1'], 
                      k_ET2=rate_constants['k_ET2'], 
                      k_EPT1=rate_constants['k_EPT1'], 
                      k_EPT2=rate_constants['k_EPT2'], 
                      reaction=reaction,
                      kinetics=kinetics,
                      )


    def calc_apparent_rate_constant(self, pH, E=None, decomposition=False):
        """
        This method calculate the apparent rate constant as a function of pH or potential. 

        Input description:
        pH (float or np.ndarray): the pH of interest
        E (float or np.ndarray): the electrode potential of interest, default value = None
        !!! NOTE: pH and E cannot be both np.ndarray
        decomposition (bool): whether to return decomposed rate constant of each channel. 

        returns:
        if decomposition == False: the overall apparent rate constant
        else: the overall and decomposed rate constants. The orders are: total, ET1, ET2, EPT1, EPT2. 
        """
        if isinstance(pH,np.ndarray) and isinstance(E,np.ndarray) and len(pH) != 1 and len(E) !=1:
            raise TypeError("pH and E cannot be both arrays")

        cH = 10**(-pH)
        chi_OxH, chi_Ox, chi_OxH3O, chi_RedH, chi_Red, chi_RedHOH = self.calc_mole_fraction(pH, E, species='all')

        if self.kinetics == 'homogeneous':
            if self.reaction == 'reduction':
                k_ET1_pH_E = chi_OxH * self.rate_constants['k_ET1']
                k_ET2_pH_E = chi_Ox * self.rate_constants['k_ET2']
                k_EPT1_pH_E = chi_OxH3O * self.rate_constants['k_EPT1']
                k_EPT2_pH_E = chi_Ox * self.rate_constants['k_EPT2']
            elif self.reaction == 'oxidation':
                k_ET1_pH_E = chi_RedH * self.rate_constants['k_ET1']
                k_ET2_pH_E = chi_Red * self.rate_constants['k_ET2']
                k_EPT1_pH_E = chi_RedH * self.rate_constants['k_EPT1']
                k_EPT2_pH_E = chi_RedHOH * self.rate_constants['k_EPT2']
            else:
                pass
        elif self.kinetics == 'Butler-Volmer' or self.kinetics == 'BV' or self.kinetics == 'Marcus-Gerischer' or self.kinetics == 'MG':
            if not(isinstance(E,float) or isinstance(E,np.ndarray) or isinstance(E,int)):
                raise ValueError("Please specify tht electrode potential. ")
            else:
                if self.reaction == 'reduction':
                    k_ET1_pH_E = chi_OxH * self.rate_constants['k_ET1'](E)
                    k_ET2_pH_E = chi_Ox * self.rate_constants['k_ET2'](E)
                    k_EPT1_pH_E = chi_OxH3O * self.rate_constants['k_EPT1'](E)
                    k_EPT2_pH_E = chi_Ox * self.rate_constants['k_EPT2'](E)
                elif self.reaction == 'oxidation':
                    k_ET1_pH_E = chi_RedH * self.rate_constants['k_ET1'](E)
                    k_ET2_pH_E = chi_Red * self.rate_constants['k_ET2'](E)
                    k_EPT1_pH_E = chi_RedH * self.rate_constants['k_EPT1'](E)
                    k_EPT2_pH_E = chi_RedHOH * self.rate_constants['k_EPT2'](E)
                else:
                    pass
        else:
            pass

        k_tot = k_ET1_pH_E + k_ET2_pH_E + k_EPT1_pH_E + k_EPT2_pH_E

        if decomposition == False:
            return k_tot
        else:
            return k_tot, k_ET1_pH_E, k_ET2_pH_E, k_EPT1_pH_E, k_EPT2_pH_E


    def calc_apparent_standard_rate_constant(self, pH):
        if self.kinetics == 'homogeneous':
            raise ValueError("The apparent standard rate constant is not defined for homogeneous systems. ")

        Eapp = self.calc_apparent_potential(pH)
        cH = 10**(-pH)

        if self.potential_dependence == False:
            KaOx = self.parameters['KaOx']
            KaRed = self.parameters['KaRed']
            K1 = self.parameters['K1']
            K2 = self.parameters['K2']
            Kw = self.parameters['Kw']
        else:
            KaOx = self.parameters['KaOx'](Eapp)
            KaRed = self.parameters['KaRed'](Eapp)
            K1 = self.parameters['K1'](Eapp)
            K2 = self.parameters['K2'](Eapp)
            Kw = self.parameters['Kw'](Eapp)

        chi_OxH = 1/(1 + KaOx/cH + K1*KaOx) 
        chi_Ox = 1/(1 + cH/KaOx + K1*cH)
        chi_OxH3O = 1/(1 + 1/KaOx/K1 + 1/K1/cH)
        chi_RedH = 1/(1 + KaRed/cH + K2*Kw/cH)
        chi_Red = 1/(1 + cH/KaRed + K2*Kw/KaRed)
        chi_RedHOH = 1/(1 + cH/K2/Kw + KaRed/K2/Kw)

        if self.reaction == 'reduction':
            k_ET1_pH_Eapp = chi_OxH * self.rate_constants['k_ET1'](Eapp)
            k_ET2_pH_Eapp = chi_Ox * self.rate_constants['k_ET2'](Eapp)
            k_EPT1_pH_Eapp = chi_OxH3O * self.rate_constants['k_EPT1'](Eapp)
            k_EPT2_pH_Eapp = chi_Ox * self.rate_constants['k_EPT2'](Eapp)
        elif self.reaction == 'oxidation':
            k_ET1_pH_Eapp = chi_RedH * self.rate_constants['k_ET1'](Eapp)
            k_ET2_pH_Eapp = chi_Red * self.rate_constants['k_ET2'](Eapp)
            k_EPT1_pH_Eapp = chi_RedH * self.rate_constants['k_EPT1'](Eapp)
            k_EPT2_pH_Eapp = chi_RedHOH * self.rate_constants['k_EPT2'](Eapp)
        else:
            pass

        ksapp = k_ET1_pH_Eapp + k_ET2_pH_Eapp + k_EPT1_pH_Eapp + k_EPT2_pH_Eapp

        return ksapp



def isnumber(dat):
    return (isinstance(dat, int) or isinstance(dat, float))

def isarray(dat):
    return (isinstance(dat, tuple) or isinstance(dat, list) or isinstance(dat, np.ndarray))

