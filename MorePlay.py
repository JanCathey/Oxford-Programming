# Code ot Calulate and Plot Radial Distribution Function (RDF) of the 3s, 3p, and 3d Orbitals

def radial_distribution_function():
    import numpy as np
    import matplotlib.pyplot as plt

    r_values = np.linspace(0, 40, 100)
    # Generates a sample RDF using a Gaussian function

    def r3s(r):
        return (r**2) * np.exp(-r / 3)

    def r3p(r):
        return (r**2) * (r**2) * np.exp(-r / 3)

    def r3d(r):
        return (r**2) * (r**4) * np.exp(-r / 3)

    rdf_3s = r3s(r_values)
    rdf_3p = r3p(r_values)
    rdf_3d = r3d(r_values)

    plt.plot(r_values, rdf_3s, label='3s Orbital', color='blue')
    plt.plot(r_values, rdf_3p, label='3p Orbital', color='orange')
    plt.plot(r_values, rdf_3d, label='3d Orbital', color='green')

    plt.title('Radial Distribution Functions of 3s, 3p, and 3d Orbitals')
    plt.xlabel('Radius (a.u.)')
    plt.ylabel('Radial Distribution Function')
    plt.legend()
    plt.grid()
    plt.show()