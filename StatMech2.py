# Performing a Computational Atomistic Thermodynamic Analysis of the adsorption behaviour of CO on Cu(III) surface using Effectie Medium Theory.
from copy import deepcopy
from ase.build import add_adsorbate, fcc111, molecule
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton

# Get energy of isolated CO molecule
CO_molec = molecule('CO')
CO_molec.calc = EMT()
Energy_CO_gas = CO_molec.get_potential_energy()

# Create a Cu(111) slab and get its energy
Cu_slab = fcc111('Cu', size=(4, 4, 2), vacuum= 10.0)
Cu_slab.calc = EMT()
energy_slab = Cu_slab.get_potential_energy()

# Run geometry optimisation of CO on Cu(III) slab and print adsorption energy
CO_ads=deepcopy(Cu_slab)
add_adsorbate(slab=CO_ads, adsorbate=CO_molec, height=3.0, position=(3.82, 2.21))
constraint = FixAtoms(mask=[atom.symbol == 'Cu' for atom in CO_ads])
CO_ads.set_constraint(constraint)
dyn = QuasiNewton(CO_ads, trajectory='CO_on_Cu111.traj')
dyn.run(fmax=0.05)
energy_CO_ads = CO_ads.get_potential_energy()

# Part A: Calculate the adsorption energy of CO on Cu(111) surface

# First, we will visualise the geometry optimisation and check the optimised structure.
from ase.visualize import view
CO_ads.set_constraint(constraint)
dyn = QuasiNewton(CO_ads, trajectory='CO_on_Cu111.traj')
dyn.run(fmax=0.05)
# view (CO_ads)

# Now, we will calculate the adsorption energy using the formula:
# E_ads,CO = E_CO+slab - (E_slab + E_CO_gas)

E_ads_CO = energy_CO_ads - (energy_slab + Energy_CO_gas)
print(f'Adsorption energy of CO on Cu(111) surface: {E_ads_CO} eV')

# Part B: Calculating the Gibbs Free Energy of Adsorption of CO on Cu(111) at 300 K and 1 bar:

import numpy as np
from ase.thermochemistry import IdealGasThermo, HarmonicThermo, HinderedThermo

vib_energies_CO_gas = np.array([0.2634])  # in eV
vib_energies_CO_ads = np.array([0.2404, 0.0827, 0.0601, 0.0600, 0.0072, 0.0065]) # in eV

thermo_CO_gas = IdealGasThermo(
    vib_energies=vib_energies_CO_gas,
    geometry='linear',
    potentialenergy=Energy_CO_gas,
    atoms=CO_molec,
    symmetrynumber=1,
    spin=0)

thermo_CO_ads = HarmonicThermo(
    vib_energies=vib_energies_CO_ads,
    potentialenergy=energy_CO_ads)

T = 300  # temperature in K
p = 1.0e+5  # pressure in Pa

g_CO_gas = thermo_CO_gas.get_gibbs_energy(temperature=T, pressure=p, verbose=False)
g_CO_ads = thermo_CO_ads.get_helmholtz_energy(temperature=T, verbose=False)
g_slab =  energy_slab  # No vibrational contribution for the slab

Pa_to_bar = 1.0e-5
adsorption_free_energy_CO = g_CO_ads - (g_slab +g_CO_gas)
print(f'Adsorption free energy of CO on Cu(11) at {T} K and {p*Pa_to_bar} bar: {adsorption_free_energy_CO}')
