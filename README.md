# PCET-pH-kinetic-model
The pH dependence of proton-coupled electron transfer (PCET) reactions is a powerful probe for elucidating their fundamental mechanisms. The PCET-pH-kinetic-model codes a general, multi-channel kinetic model for describing the pH dependence of PCET reactions that can be applied to both homogeneous and electrochemical reactions.  

## Installation 
To use the kinetic model, simply download the code and add it to your `$PYTHONPATH` variable.

## Documentation

### Initialization

#### General consideration
To start a calculation, create a `pHKineticModel` object and input the parameters: 
```python
from PCET_pH import pHKineticModel

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
```
The `pHKineticModel` object stores all thermodyanmics and kinetics paramaters in the PCET kinetic model, which will be used to calculate the mole fractions, apparent potentials and apparent rate constants as functions of pH (and electrode potential). 

The requiered parameters are listed below:

1. `KaOx` (float or function): acid dissociation equilibrium constant of the oxidized species.

2. `KaRed` (float or function): acid dissociation equilibrium constant of the reduced species.

3. `K1` (float or function): association equilibrium constant of the oxidized species with H3O+.

4. `K2` (float or function): association equilibrium constant of the reduced species with OH-.

5. `E_ET1` (float or function): redox potential of the ET1 channel (in Volts).
    
6. `k_ET1`/`k_ET2`/`k_EPT1`/`k_EPT2` (float or tuple): rate constants (and other parameters) of the ET1/ET2/EPT1/EPT2 channel
This class can be used to describe both homogeneous and electrochemical reactions. For the later, only Butler-Volmer or Marcus-Gerischer kinetic models are supported. The input format of these `k` parameters varies with the kinetic model choice. 
```
  if reaction = 'homogeneous': k = rate constant (float)

  if reaction = 'Butler-Volme': k = (standard rate constant, cathodic charge transfer coefficient)

  if reaction = 'Marcus-Gerischer': k = (maximum rate constant, reorganization energy (in eV))
```                                               
7. `kinetics` (string): choose among `'homogeneous'`, `'Butler-Volmer'` or `'BV'`, and `'Marcus-Gerischer'` or `'MG'`

8. `reaction` (string): choose between `'reduction'` and `'oxidation'`

These parameters can be summarized in the generalized square diagram of PCET. Detailed definitions can be found in Ref. 1. 
<div align="center">
  <img src="./PCET_diagram.jpg" alt="PCET" width="600">
</div>

#### Potential-dependent equilibrium constants
This class also enables the use of potential-dependent equilibrium constants for electrochemical systems. This can be done by setting the `K` parameters to be user-defined functions. For example: 

```python
# define the potential-dependent thermodynamics parameters
def E_ET1_E(E):
    return 0.12*(E-1.1)+1.1

def KaOx_E(E):
    return 10**(-(-0.12*(E-1.1)/0.059+1.2))

# create an instance of the kinetic model and use the pre-defined functions as input parameters
system = pHKineticModel(KaOx=KaOx_E, 
                        KaRed=1e-9, 
                        K1=1e1, 
                        K2=1e3, 
                        E_ET1=E_ET1_E, 
                        k_ET1=(36,0.37), 
                        k_ET2=(12,0.37), 
                        k_EPT1=(12,1.01), 
                        k_EPT2=(12,1.01), 
                        reaction='reduction', 
                        kinetics='MG',
                        )
```


### Examples
