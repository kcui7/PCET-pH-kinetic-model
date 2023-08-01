import numpy as np
import matplotlib.pyplot as plt
from PCET_pH import pHKineticModel

# create an instance of the kinetic model and input parameters
system = pHKineticModel(KaOx=1e-3, 
                        KaRed=1e-11, 
                        K1=1e2, 
                        K2=1e4, 
                        E_ET1=1.0, 
                        k_ET1=(1e6,0.5), 
                        k_ET2=(1e4,0.5), 
                        k_EPT1=(1e5,0.5), 
                        k_EPT2=(1e7,0.5), 
                        reaction='reduction', 
                        kinetics='BV',
                        )


#===============================================================================
# Example 3.1. calculate apparent standard rate constants as functions of pH
#===============================================================================

# define the pH of interest
pH = np.linspace(0,14,100)

# calculate the apparent standard rate constant by calling
# system.calc_apparent_standard_rate_constant(pH)
plt.semilogy(pH, system.calc_apparent_standard_rate_constant(pH), 'r-', lw=1.5, label=r'All $\alpha=0.5$')

# reset charge transfer coefficients
system.set_parameters(k_ET1=(1e6,0.4), 
                      k_ET2=(1e4,0.4), 
                      k_EPT1=(1e5,0.4), 
                      k_EPT2=(1e7,0.4), 
                      )
plt.semilogy(pH, system.calc_apparent_standard_rate_constant(pH), 'b-', lw=1.5, label=r'All $\alpha=0.4$')

# reset charge transfer coefficients
system.set_parameters(k_ET1=(1e6,0.6), 
                      k_ET2=(1e4,0.6), 
                      k_EPT1=(1e5,0.6), 
                      k_EPT2=(1e7,0.6), 
                      )
plt.semilogy(pH, system.calc_apparent_standard_rate_constant(pH), '-', c=(1,0,1), lw=1.5, label=r'All $\alpha=0.6$')

# parameters used to make a fancier plot
plt.legend(fontsize=15)
plt.xlabel('pH', fontsize=17)
plt.ylabel(r'$k_{\rm s,app}\ /\ \rm s^{-1}$', fontsize=17)
plt.xlim(0,14)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.show()
