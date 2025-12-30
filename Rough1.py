'''
Rough Code for Radial Distribution Function Calculation
First, we will Calculate RDFs for 1s, 2s, and 2p Orbitals, then normalise them so the maximum value of each RDF is 1.
Note: The RDF is 4pi*r^2*|R(r)|^2 where R(r) is the radial wavefunction.
'''

import numpy as np
import matplotlib.pyplot as plt

r_values = np.linspace(0, 20, 200)

# Define Unnormalised Wavefunctions
def psi_1s(r):
    psi_1s = np.exp(-r/2)
    return psi_1s

def psi_2s(r):
    psi_2s = (2 - r)*np.exp(-r/2) 
    return psi_2s

def psi_2p(r):
    psi_2p = (6 - 6*r + r**2) * np.exp(-r/2)
    return psi_2p

# Calculate Radial Distribution Functions
def rdf(psi, r):
   return 4 * np.pi * r**2 * np.abs(psi(r))**2

rdf_1s = rdf(psi_1s, r_values)
rdf_2s = rdf(psi_2s, r_values)
rdf_2p = rdf(psi_2p, r_values)

# Normalise Radial Distribution Functions
def normalise_rdf(rdf):
    max_value = np.max(rdf)
    return rdf / max_value 

rdf_1s_normalised = normalise_rdf(rdf_1s)
rdf_2s_normalised = normalise_rdf(rdf_2s)
rdf_2p_normalised = normalise_rdf(rdf_2p)

# Plotting the Normalised RDFs
plt.figure(figsize=(10, 6))
plt.plot(r_values, rdf_1s_normalised, label='1s Orbital', color='blue')
plt.plot(r_values, rdf_2s_normalised, label='2s Orbital', color='orange')
plt.plot(r_values, rdf_2p_normalised, label='2p Orbital', color='green')
plt.title('Normalised Radial Distribution Functions')
plt.xlabel('Radius (r)')
plt.ylabel('Normalised RDF')
plt.legend()
plt.grid()
plt.show()  

# Now, we will write a routine to Numerically Integrate the RDFs.
# We will use the Trapezoidal Rule for Numerical Integration.
def integrate_rdf(rdf, r):
    integral = np.trapezoid(rdf, r)
    return integral

integral_1s = integrate_rdf(rdf_1s, r_values)
integral_2s = integrate_rdf(rdf_2s, r_values)
integral_2p = integrate_rdf(rdf_2p, r_values)
print(f'Integral of RDF for 1s Orbital: {integral_1s}')
print(f'Integral of RDF for 2s Orbital: {integral_2s}')
print(f'Integral of RDF for 2p Orbital: {integral_2p}')
