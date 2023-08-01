import numpy as np
import matplotlib.pyplot as plt
from PCET_pH import pHKineticModel

# create an instance of the kinetic model and input parameters
system = pHKineticModel(KaOx=1e-3, 
                        KaRed=1e-11, 
                        K1=1e2, 
                        K2=1e4, 
                        E_ET1=1.0, 
                        k_ET1=1e4, 
                        k_ET2=1e5, 
                        k_EPT1=1e6, 
                        k_EPT2=1e7, 
                        reaction='reduction', 
                        kinetics='homogeneous',
                        )


#===============================================================================
# Example 1.1. calculate mole fractions as functions of pH
#===============================================================================

# define the pH of interest
pH = np.linspace(0,14,100)

# calculate the mole fraction by calling
# system.calc_mole_fraction(pH)
# By default this method returns all mole fractions in the order of 
# OxH, Ox, OxH3O, RedH, Red, RedHOH
chi_OxH, chi_Ox, chi_OxH3O, chi_RedH, chi_Red, chi_RedHOH = system.calc_mole_fraction(pH)

# plot mole fractions of the oxidized species
plt.plot(pH, chi_OxH, 'k-', lw=1.5, label=r'$\chi_{\rm OxH^+}$')
plt.plot(pH, chi_Ox, 'b-', lw=1.5, label=r'$\chi_{\rm Ox}$')
plt.plot(pH, chi_OxH3O, 'r-', lw=1.5, label=r'$\chi_{\rm Ox\cdot\!\cdot\!\cdot H_3O^+}$')

# parameters used to make a fancier plot
plt.xlim(0,14)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$\chi_{{\rm Ox}_i}$', fontsize=17)
plt.legend(loc=(0.63,0.65),fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()


# plot mole fractions of the oxidized species
plt.plot(pH, chi_RedH, 'k-', lw=1.5, label=r'$\chi_{\rm RedH}$')
plt.plot(pH, chi_Red, 'b-', lw=1.5, label=r'$\chi_{\rm Red^-}$')
plt.plot(pH, chi_RedHOH, '-', lw=1.5, c=(1,0,1), label=r'$\chi_{\rm RedH\cdot\!\cdot\!\cdot OH^-}$')

# parameters used to make a fancier plot
plt.xlim(0,14)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$\chi_{{\rm Red}_i}$', fontsize=17)
plt.legend(loc=(0.03,0.65), fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#===============================================================================
# Example 1.2. calculate apparent potential as functions of pH
#===============================================================================

# define the pH of interest
pH = np.linspace(0,14,100)

# calculate the apparent potential by calling
# system.calc_apparent_potential(pH)
plt.plot(pH, system.calc_apparent_potential(pH), 'r-', lw=1.5, label=r'$E^{\rm o}_{\rm app}$')

# parameters used to make a fancier plot
plt.xlim(0,14)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$E$ / V vs. NHE', fontsize=17)
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#===============================================================================
# Example 1.3. calculate apparent rate constants as functions of pH
#===============================================================================

# define the pH of interest
pH = np.linspace(0,14,100)

# calculate the apparent rate constant for reduction by calling
# system.calc_apparent_potential(pH)
plt.semilogy(pH, system.calc_apparent_rate_constant(pH), 'r-', lw=1.5, label=r'$k^{\rm Red}_{\rm app}$')

# calculate the apparent rate constant for oxidation
# assume the redox potential of the external molecular reactant is EM = 0.85 V

EM = 0.85
F = 96485
R = 8.314
T = 298

# calculate the oxidation rate constants of each channel using detailed balance relation
# k_red/k_ox = ln( -F(EM-E)/RT )
k_ET1_ox = np.exp(F*(EM-system.get_parameters('E_ET1'))/R/T) * system.get_parameters('k_ET1')
k_ET2_ox = np.exp(F*(EM-system.get_parameters('E_ET2'))/R/T) * system.get_parameters('k_ET2')
k_EPT1_ox = np.exp(F*(EM-system.get_parameters('E_EPT1'))/R/T) * system.get_parameters('k_EPT1')
k_EPT2_ox = np.exp(F*(EM-system.get_parameters('E_EPT2'))/R/T) * system.get_parameters('k_EPT2')

# update the parameters in the kinetic model
system.set_parameters(k_ET1=k_ET1_ox,
                      k_ET2=k_ET2_ox,
                      k_EPT1=k_EPT1_ox,
                      k_EPT2=k_EPT2_ox,
                      reaction='oxidation',
                      )

# calculate apparent rate constant
plt.semilogy(pH, system.calc_apparent_rate_constant(pH), 'b-', lw=1.5, label=r'$k^{\rm Ox}_{\rm app}$')

# parameters used to make a fancier plot
plt.xlim(0,14)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$k_{\rm app}$ / $\rm M^{-1}\,s^{-1}$', fontsize=17)
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()
