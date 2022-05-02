import sys, os
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
from sys import stdout
import numpy as np
import mdtraj
import grand

# Load PDB
pdb = PDBFile('3rlp-uvt1.pdb')

# Create system
forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=10.0*angstrom,
                                 switchDistance=8.0*angstrom, constraints=HBonds)


# Define the barostat
system.addForce(MonteCarloBarostat(1*bar, 286*kelvin, 25))

# Define integrator
integrator = BAOABIntegrator(286*kelvin, 1.0/picosecond, 0.002*picosecond)

# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object
simulation = Simulation(pdb.topology, system, integrator, platform)

# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(286*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Add a reporter - write out every 5 ps
simulation.reporters.append(StateDataReporter(stdout, 2500, step=True, time=True, potentialEnergy=True,
                                              temperature=True, volume=True))

simulation.reporters.append(DCDReporter('3rlp-npt.dcd', 2500)) 
# Run for 500 ps
simulation.step(250000)

# Write out a PDB
positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
pdb.topology.setPeriodicBoxVectors(simulation.context.getState(getPositions=True).getPeriodicBoxVectors())
with open('3rlp-npt.pdb', 'w') as f:
    PDBFile.writeFile(pdb.topology, positions, f)




