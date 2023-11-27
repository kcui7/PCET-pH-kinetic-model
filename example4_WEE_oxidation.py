import numpy as np
import matplotlib.pyplot as plt
from PCET_pH_Trp import pHKineticModelTrp

# create an instance of the kinetic model and input parameters
system = pHKineticModelTrp(KaNH3=10**(-7.5), 
                           KOH_NH3=7e4, 
                           KOH_NH2=4.2e2, 
                           k_ET1_NH3=7.1e5/2, 
                           k_EPT1_NH3=7.1e5/2, 
                           k_EPT2_NH3=4e9, 
                           k_ET1_NH2=7.1e5/2, 
                           k_EPT1_NH2=7.1e5/2, 
                           k_EPT2_NH2=4e8, 
                           )

# define the pH of interest
pH = np.linspace(0,14,100)

# use a different set of parameters
#system.set_parameters(KOH_NH3=7e5,
#                      KOH_NH2=4.2e3,
#                      k_EPT2_NH3=4e8)


#===============================================================================
# Example 4.1. calculate mole fractions as functions of pH
#===============================================================================

# calculate the mole fraction by calling
# system.calc_mole_fraction(pH)
# By default this method returns all mole fractions in the order of 
# TrpNH_NH3, TrpNH_NH2, TrpNH_NH3_OH, TrpNH_NH3_OH
chi_TrpNH_NH3, chi_TrpNH_NH2, chi_TrpNH_NH3_OH, chi_TrpNH_NH2_OH = system.calc_mole_fraction(pH, species='all')

# plot mole fractions of the oxidized species
plt.plot(pH, chi_TrpNH_NH3, 'k-', lw=1.5, label=r'$\chi_{\rm Trp(NH_3^+\!)NH}$')
plt.plot(pH, chi_TrpNH_NH2, 'b--', lw=1.5, label=r'$\chi_{\rm Trp(NH_2)NH}$')
plt.plot(pH, chi_TrpNH_NH3_OH, 'g-', lw=1.5, label=r'$\chi_{\rm Trp(NH_3^+\!)NH\cdot\!\cdot\!\cdot OH^-}$')
plt.plot(pH, chi_TrpNH_NH2_OH, 'r--', lw=1.5, label=r'$\chi_{\rm Trp(NH_2)NH\cdot\!\cdot\!\cdot OH^-}$')

# parameters used to make a fancier plot
plt.xlim(0,14)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$\chi_{{\rm Red}_i}$', fontsize=17)

plt.legend(loc=(0.01,0.35), fontsize=16, frameon=False)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#===============================================================================
# Example 4.2. calculate apparent rate constants as functions of pH
#===============================================================================

# Experimental data for WEE oxidation
pH_exp_WEE = np.array([2.0332750242954325, 2.177010689990282, 3.068345966958212, 4.9018931000971815, 5.586487852283771, 6.096746355685131, 6.447409135082604, 6.680707482993197, 7.0472303206997085, 7.602829931972789, 8.646421768707484, 10.518903790087464, 11.36806997084548])
logk_exp_WEE = np.array([5.866705580932532, 5.856030274776504, 5.901650069340521, 5.803907357930283, 6.288213678495037, 6.611587310292954, 6.841625467167067, 6.98110960503326, 7.295629588422876, 7.7925112237622525, 8.195567691234888, 8.436222391111878, 8.539508426636162])

# calculate the apparent rate constant for oxidation by calling
# system.calc_apparent_potential(pH)
plt.semilogy(pH, system.calc_apparent_rate_constant(pH), '-', lw=1.5, label=r'$k_{\rm app}^{\rm Ox}$', c='#F981FF')

# plot experimental data
plt.semilogy(pH_exp_WEE[:-4], 10**(logk_exp_WEE)[:-4], '^', ms=7, c='#F981FF')
plt.semilogy(pH_exp_WEE[-4:], 10**(logk_exp_WEE)[-4:], 'd', ms=7, c='#FE9401')

plt.xlim(0,14)
plt.ylim(2e5,2e9)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$k_{\rm app}$ /$\rm M^{-1}\,s^{-1}$', fontsize=17)

plt.legend(loc=2,fontsize=15,framealpha=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()



#===============================================================================
# Example 4.3. calculate contribution of each channel to kapp as functions of pH
#===============================================================================

# calculate the desomposed apparent rate constant for oxidation by calling
# system.calc_apparent_potential(pH, decomposition=True)
k_tot, k_ET1_NH3, k_EPT1_NH3, k_EPT2_NH3, k_ET1_NH2, k_EPT1_NH2, k_EPT2_NH2 = system.calc_apparent_rate_constant(pH, decomposition=True)

plt.plot(pH, (k_ET1_NH3 + k_EPT1_NH3)/k_tot, '-', lw=1.5, label=r'ET1+EPT1', alpha=1, c='k')
plt.plot(pH, (k_ET1_NH2 + k_EPT1_NH2)/k_tot, '--', lw=1.5, label=r"ET1'+EPT1'", alpha=1, c='b')
plt.plot(pH, k_EPT2_NH3/k_tot, '-', lw=1.5, label=r'EPT2', alpha=1, c='g')
plt.plot(pH, k_EPT2_NH2/k_tot, '--', lw=1.5, label=r"EPT2'", alpha=1, c='r')

plt.xlim(0,14)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'Contribution to $k_{\rm app}$', fontsize=17)

plt.legend(fontsize=14.5,frameon=False)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()

