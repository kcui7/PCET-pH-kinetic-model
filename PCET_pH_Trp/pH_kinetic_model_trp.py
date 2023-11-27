import numpy as np

class GeneralizedSquareDiagramTrp(object):
    """
    This class stores all thermodynamics parameters which will be used to calculate the mole fractions at given pH. 
    This class is an modification of the original GeneralizedSquareDiagram class in the following repository:
        https://github.com/kcui7/PCET-pH-kinetic-model.git
    This modification is tailored to model the pH dependence in Trp-derivative systems. 
    Due to the complexity, the potential-dependence feature is disabled. 
    In addition, we focus on the OXIDATION reaction for this system, so the mole fractions of oxidized species and the apparent potential will NOT be calculated

    Input description and data type:
    KaNH3 (float): acid dissociation equilibrium constant of the NH3+ group 
    KOH_NH3 (float): association equilibrium constant of the reduced species Trp(NH3+)NH with OH-
    KOH_NH2 (float): association equilibrium constant of the reduced species Trp(NH2)NH with OH-
    T (float): the temperature, default value = 298K
    Kw (float): water self-dissociation equilibrium constant at the given temperature, default value = 1e-14
    """
    def __init__(self, KaNH3=1e-7, KOH_NH3=7e4, KOH_NH2=4e3, T=298, Kw=1e-14):
        self.T = T
        self.R = 8.314
        self.F = 96485
        self.parameters = {'KaNH3':KaNH3,
                           'KOH_NH3':KOH_NH3,
                           'KOH_NH2':KOH_NH2,
                           'Kw':Kw,
                           }

    def __repr__(self):
        str = "The following parameters are defined in the generlized square diagram\n"
        for parameter in ['KaNH3', 'KOH_NH3', 'KOH_NH2', 'Kw']:
            str += "%s = %.2e\n"%(parameter, self.parameters[parameter])

        return str


    def get_parameters(self, key='all'):
        if key == 'all':
            return self.parameters
        else:
            return self.parameters[key]


    def set_parameters(self, **kwargs): 
        # re-initialize the instance
        for key in kwargs.keys():
            if key in self.parameters.keys() and key.startswith('K'):
                self.parameters[key] = kwargs[key]
            elif key == 'T':
                self.T = kwargs[key]
            else:
                pass


    def calc_mole_fraction(self, pH, species='all'):
        """
        This method calculate the mole fraction of the specified species as a function of pH

        Input description:
        pH (float or np.ndarray): the pH of interest
        species (string): the species of interest, possible choices: 'TrpNH_NH3', 'TrpNH_NH2', 'TrpNH_NH3_OH', 'TrpNH_NH2_OH'

        returns:
        if species == 'all': mole fractions of all species, with the order of: TrpNH_NH3, TrpNH_NH2, TrpNH_NH3_OH, TrpNH_NH2_OH
        else: mole fraction of the specific species. 
        """
        cH = 10**(-pH)

        KaNH3 = self.parameters['KaNH3']
        KOH_NH3 = self.parameters['KOH_NH3']
        KOH_NH2 = self.parameters['KOH_NH2']
        Kw = self.parameters['Kw']

        chi_TrpNH_NH3 = 1/(1 + KaNH3/cH + KOH_NH3*Kw/cH + KOH_NH2*Kw*KaNH3/cH/cH)
        chi_TrpNH_NH2 = (KaNH3/cH)/(1 + KaNH3/cH + KOH_NH3*Kw/cH + KOH_NH2*Kw*KaNH3/cH/cH)
        chi_TrpNH_NH3_OH = (KOH_NH3*Kw/cH)/(1 + KaNH3/cH + KOH_NH3*Kw/cH + KOH_NH2*Kw*KaNH3/cH/cH)
        chi_TrpNH_NH2_OH = (KOH_NH2*Kw*KaNH3/cH/cH)/(1 + KaNH3/cH + KOH_NH3*Kw/cH + KOH_NH2*Kw*KaNH3/cH/cH)

        if species == 'all':
            return chi_TrpNH_NH3, chi_TrpNH_NH2, chi_TrpNH_NH3_OH, chi_TrpNH_NH2_OH
        elif species == 'TrpNH_NH3':
            return chi_TrpNH_NH3
        elif species == 'TrpNH_NH2':
            return chi_TrpNH_NH2
        elif species == 'TrpNH_NH3_OH':
            return chi_TrpNH_NH3_OH
        elif species == 'TrpNH_NH2_OH':
            return chi_TrpNH_NH2_OH
        else:
            raise KeyError("'%s'"%species)




class pHKineticModelTrp(GeneralizedSquareDiagramTrp):
    """
    This class stores all thermodynamics and kinetics parameters in the PCET kinetic model which can be used to 
    calculate the mole fractions and apparent rate constants as functions of pH. 
    This class is an modification of the original pHKineticModel class in the following repository:
        https://github.com/kcui7/PCET-pH-kinetic-model.git
    This modification is tailored to model the pH dependence in Trp-derivative systems. 
    This class can only be used to describe both HOMOGENEOUS reactions. 

    Input description and data type:
    KaNH3 (float): acid dissociation equilibrium constant of the NH3+ group 
    KOH_NH3 (float): association equilibrium constant of the reduced species Trp(NH3+)NH with OH-
    KOH_NH2 (float): association equilibrium constant of the reduced species Trp(NH2)NH with OH-
    T (float): the temperature, default value = 298K
    Kw (float): water self-dissociation equilibrium constant at the given temperature, default value = 1e-14
    
    k_ET1_NH3/k_EPT1_NH3/k_EPT2_NH3/k_ET1_NH2/k_EPT1_NH2/k_EPT2_NH2 (float): rate constants of the ET1/EPT1/EPT2/ET1'/EPT1'/EPT2' channels
    """
    def __init__(self, KaNH3=1e-7, KOH_NH3=7e4, KOH_NH2=4e3, T=298, Kw=1e-14,
                       k_ET1_NH3=1e5, k_EPT1_NH3=1e5, k_EPT2_NH3=1e5,
                       k_ET1_NH2=1e5, k_EPT1_NH2=1e5, k_EPT2_NH2=1e5):
        
        super().__init__(KaNH3=KaNH3, KOH_NH3=KOH_NH3, KOH_NH2=KOH_NH2, T=T, Kw=Kw)

        self.rate_constants = {}

        self.rate_constants['k_ET1_NH3'] = k_ET1_NH3
        self.rate_constants['k_EPT1_NH3'] = k_EPT1_NH3
        self.rate_constants['k_EPT2_NH3'] = k_EPT2_NH3
        self.rate_constants['k_ET1_NH2'] = k_ET1_NH2
        self.rate_constants['k_EPT1_NH2'] = k_EPT1_NH2
        self.rate_constants['k_EPT2_NH2'] = k_EPT2_NH2


    def __repr__(self):
        str = super().__repr__()

        str2 = "\nThe following rate constants are defined for the homogeneous system. \n"
        for key in ['k_ET1_NH3', 'k_EPT1_NH3', 'k_EPT2_NH3', 'k_ET1_NH2', 'k_EPT1_NH2','k_EPT2_NH2']:
            str2 += "%s: %e\n"%(key, self.rate_constants[key])
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
        if key.startswith('K'):
            return self.parameters[key]
        elif key.startswith('k'):
            return self.rate_constants[key]
        else:
            raise KeyError("'%s'"%key)


    def set_parameters(self, **kwargs):
        # re-initialize the instance
        for key in kwargs.keys():
            if key in self.parameters.keys() and key.startswith('K'):
                self.parameters[key] = kwargs[key]
            elif key == 'T':
                self.T = kwargs[key]
            elif key in self.rate_constants.keys():
                self.rate_constants[key] = kwargs[key]
            else:
                pass


    def calc_apparent_rate_constant(self, pH, decomposition=False):
        """
        This method calculate the apparent rate constant as a function of pH.

        Input description:
        pH (float or np.ndarray): the pH of interest
        decomposition (bool): whether to return decomposed rate constant of each channel. 

        returns:
        if decomposition == False: the overall apparent rate constant
        else: the overall and decomposed rate constants. The orders are: total, ET1, EPT1, EPT2, ET1', EPT1', EPT2'. 
        """

        cH = 10**(-pH)
        chi_TrpNH_NH3, chi_TrpNH_NH2, chi_TrpNH_NH3_OH, chi_TrpNH_NH2_OH = self.calc_mole_fraction(pH, species='all')

        k_ET1_NH3_pH = chi_TrpNH_NH3 * self.rate_constants['k_ET1_NH3']
        k_EPT1_NH3_pH = chi_TrpNH_NH3 * self.rate_constants['k_EPT1_NH3']
        k_EPT2_NH3_pH = chi_TrpNH_NH3_OH * self.rate_constants['k_EPT2_NH3']
        k_ET1_NH2_pH = chi_TrpNH_NH2 * self.rate_constants['k_ET1_NH2']
        k_EPT1_NH2_pH = chi_TrpNH_NH2 * self.rate_constants['k_EPT1_NH2']
        k_EPT2_NH2_pH = chi_TrpNH_NH2_OH * self.rate_constants['k_EPT2_NH2']

        k_tot = k_ET1_NH3_pH + k_EPT1_NH3_pH + k_EPT2_NH3_pH + k_ET1_NH2_pH + k_EPT1_NH2_pH + k_EPT2_NH2_pH

        if decomposition == False:
            return k_tot
        else:
            return k_tot, k_ET1_NH3_pH, k_EPT1_NH3_pH, k_EPT2_NH3_pH, k_ET1_NH2_pH, k_EPT1_NH2_pH, k_EPT2_NH2_pH

