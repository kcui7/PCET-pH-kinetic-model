import numpy as np
import matplotlib.pyplot as plt
from PCET_pH import pHKineticModel

# define the potential-dependent thermodynamics parameters
def E_ET1_E(E):
    return 0.12*(E-1.1)+1.1

def KaOx_E(E):
    return 10**(-(-0.12*(E-1.1)/0.059+1.2))

def KaRed_E(E):
    return 10**(-(-0.12*(E-1.1)/0.059+8.2))


# create an instance of the kinetic model and input parameters
system = pHKineticModel(KaOx=KaOx_E, 
                        KaRed=KaRed_E, 
                        K1=10**(0.2), 
                        K2=1e3, 
                        E_ET1=E_ET1_E, 
                        k_ET1=(36,0.37), 
                        k_ET2=(12,0.37), 
                        k_EPT1=(12,1.01), 
                        k_EPT2=(12,1.01), 
                        reaction='reduction', 
                        kinetics='MG',
                        )


#===============================================================================
# Example 2.1. calculate apparent rate constants as functions of applied potential
#===============================================================================

# plot experimental data
E_exp_pH1 = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
k_exp_pH1 = [36,36,34,31,28,24,19,12,6,2]
E_exp_pH3 = [-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
k_exp_pH3 = [33,32,28,21,13,7,3,1,0.3,0.1]
E_exp_pH5 = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7]
k_exp_pH5 = [12,12,10,7,3,0.9,0.2,0.04]

plt.plot(E_exp_pH1, k_exp_pH1, 'ro', ms=7)
plt.plot(E_exp_pH3, k_exp_pH3, 'bx', ms=7)
plt.plot(E_exp_pH5, k_exp_pH5, 's', ms=7, color=(1,0,1))

# define the potential of interest
E = np.linspace(1.1,-0.5,100)

# calculate the apparent rate constant by calling
# system.calc_apparent_rate_constant(pH, E)
# this returns the total apparent rate constant at the given pH and E. 
plt.plot(E, system.calc_apparent_rate_constant(pH=1, E=E), 'r-', lw=1.5, label='pH=1')
plt.plot(E, system.calc_apparent_rate_constant(pH=3, E=E), 'b-', lw=1.5, label='pH=3')
plt.plot(E, system.calc_apparent_rate_constant(pH=5, E=E), '-', lw=1.5, label='pH=5', color=(1,0,1))

# parameters used to make a fancier plot
plt.legend(fontsize=15)
plt.xlabel(r'$E$ / V vs. NHE', fontsize=17)
plt.ylabel(r'$k_{\rm app}\ /\ \rm s^{-1}\times 10^{6}$', fontsize=17)
plt.xlim(1.1,-0.2)
plt.ylim(-1,38)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#===============================================================================
# Example 2.2. calculate apparent maximum rate constant as a function of pH
#===============================================================================

# define the pH of interest
pH = np.linspace(0,14,100)

# Einf is a very negative potential where rate constants plateau
Einf = -0.1

# calculate the apparent rate constant at Einf as a function of pH
plt.plot(pH, system.calc_apparent_rate_constant(pH=pH, E=Einf), 'k-', lw=1.5, label=r'$k_{\rm max,app}^{\rm Red}$')

# plot experimental data
pH_data = [1,3,5]
kapp_data = [36,33,12]
plt.plot(pH_data, kapp_data, 'ks', ms=7, mec='k', mfc='w', mew=1.5, label='Experiment')

# parameters used to make a fancier plot
plt.legend(loc=3,fontsize=15,framealpha=1)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$k_{\rm max, app}$ / ${\rm s}^{-1}\times 10^6$', fontsize=17)
plt.xlim(0,14)
plt.ylim(0,42)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#===============================================================================
# Example 2.3. decompose the total apparent rate constant
#===============================================================================

# plot experimental data
E_exp_pH3 = [-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
k_exp_pH3 = [33,32,28,21,13,7,3,1,0.3,0.1]
plt.plot(E_exp_pH3, k_exp_pH3, 'bx', ms=7)

# define the potential of interest
E = np.linspace(1.1,-0.5,100)

# calculate the decomposed apparent rate constant by calling
# system.calc_apparent_rate_constant(pH, E, decomposition=True)
# this returns the total and decomposed apparent rate constants at the given pH and E. 
# the order is k_tot, k_ET1, k_ET2, k_EPT1, k_EPT2
k_tot, k_ET1, k_ET2, k_EPT1, k_EPT2 = system.calc_apparent_rate_constant(pH=3, E=E, decomposition=True)

# plot the rate constants
plt.plot(E, k_tot, 'b-', lw=1.5, label='total')
plt.plot(E, k_ET1, 'k-', lw=1.5, label='ET1')
plt.plot(E, k_ET2, 'k--', lw=1.5, label='ET2')
plt.plot(E, k_EPT1, 'k-.', lw=1.5, label='EPT1')
plt.plot(E, k_EPT2, 'k:', lw=1.5, label='EPT2')

# parameters used to make a fancier plot
plt.legend(fontsize=15)
plt.xlabel(r'$E$ / V vs. NHE', fontsize=17)
plt.title(r'pH = 3', fontsize=17)
plt.ylabel(r'$k_{\rm app}\ /\ \rm s^{-1}\times 10^{6}$', fontsize=17)
plt.xlim(1.1,-0.2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#===============================================================================
# Example 2.4. calculate mole fractions as functions of E
#===============================================================================

# define the potential of interest
E = np.linspace(1.1,-0.5,100)

# calculate the mole fraction by calling 
# system.calc_mole_fraction(pH, E)
chi_OxH, chi_Ox, chi_OxH3O, chi_RedH, chi_Red, chi_RedHOH = system.calc_mole_fraction(pH=3, E=E)

# plot mole fractions of the oxidized species
plt.plot(E, chi_OxH, 'k-', lw=1.5, label=r'$\chi_{\rm OxH^+}$')
plt.plot(E, chi_Ox, 'k--', lw=1.5, label=r'$\chi_{\rm Ox}$')
plt.plot(E, chi_OxH3O, 'k-.', lw=1.5, label=r'$\chi_{\rm Ox\cdot\!\cdot\!\cdot H_3O^+}$')

# parameters used to make a fancier plot
plt.xlim(1.1,-0.2)
plt.title(r'pH = 3', fontsize=17)
plt.xlabel(r'$E$ / V vs. NHE', fontsize=17)
plt.ylabel(r'$\chi_{{\rm Ox}_i}$', fontsize=17)
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()
