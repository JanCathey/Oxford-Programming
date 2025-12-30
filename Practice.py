# Preliminary code for writing Radial Distribtion Functions (RDF)
import numpy as np
import matplotlib.pyplot as plt

r_values = np.linspace(0, 30, 100)
# Generates a sample RDF using a Gaussian function

def r3s(r):
    """
    Calculates the value of the radial wavefunction of the 3s orbital for a given radius r.
    
    Args:
        r (float): The radius at which to evaluate the wavefunction.
    
    Returns:
        
        (float): The value of the radial wavefunction at radius r.
    """
    
    rs3 = 2/(3*np.sqrt(3)) * (r/3) * (2 - r/3) * np.exp(-r/3)
    return rs3

def r3p(r):
    """
    Calculates the value of the radial wavefunction of the 3p orbital for a given radius r.
    
    Args:
        r (float): The radius at which to evaluate the wavefunction.
    
    Returns:
        (float): The value of the radial wavefunction at radius r.
    """
    
    rp3 = (2*np.sqrt(2))/(9*np.sqrt(3)) * (((2*r)/3) * (1 - (r/6))) * np.exp(-r/3)
    return rp3

def r3d(r):
    """
    Calculates the value of the radial wavefunction of the 3d orbital for a given radius r.
    
    Args:
        r (float): The radius at which to evaluate the wavefunction.
    Returns:
        (float): The value of the radial wavefunction at radius r.
    """
    
    rd3 = (4/(81*np.sqrt(30))) * ((r/3)**2) * np.exp(-r/3)
    return rd3

# Calculate R(r) from r3s, r3p, r3d functions
wavefunction_3s = r3s(r_values)
wavefunction_3p = r3p(r_values)
wavefunction_3d = r3d(r_values)

# Plotting the radial wavefunctions
plt.plot(r_values, wavefunction_3s, label='3s Orbital')
plt.plot(r_values, wavefunction_3p, label='3p Orbital')
plt.plot(r_values, wavefunction_3d, label='3d Orbital')
plt.title('Radial Wavefunctions of 3s, 3p, and 3d Orbitals')
plt.xlabel('Radius (r)')
plt.ylabel('Radial Wavefunction R(r)')
plt.legend()
plt.grid()
plt.show()  

# Now we can proceed to calculate and plot the Radial Distribution Function (RDF)
def radial_distribution_function(r, wavefunction):
    """
    Calculates the Radial Distribution Function (RDF) for a given radial wavefunction.
    
    Args:
        r (numpy.ndarray): Array of radius values.
        wavefunction (numpy.ndarray): Corresponding radial wavefunction values.
    
    Returns:
        float: Radial Distribution Function values at r.
    """
    rdf = 4 * np.pi * r**2 * wavefunction**2
    return rdf
rdf_3s = radial_distribution_function(r_values, wavefunction_3s)
rdf_3p = radial_distribution_function(r_values, wavefunction_3p)
rdf_3d = radial_distribution_function(r_values, wavefunction_3d)

# Plotting the Radial Distribution Functions
plt.plot(r_values, rdf_3s, label='RDF 3s Orbital')
plt.plot(r_values, rdf_3p, label='RDF 3p Orbital')
plt.plot(r_values, rdf_3d, label='RDF 3d Orbital')
plt.title('Radial Distribution Functions of 3s, 3p, and 3d Orbitals')
plt.xlabel('Radius (r)')
plt.ylabel('Radial Distribution Function RDF(r)')
plt.legend()
plt.grid()
plt.show()

