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
pdb = PDBFile('1c5t-ghosts.pdb')

# Load restart file
rst7 = AmberInpcrdFile('1c5t-prod.rst7')

# Shouldn't need to add ghosts as these can just be read in from before (all frames contained below)
ghosts = grand.utils.read_ghosts_from_file('1c5t-prod.txt')

ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=10.0*angstroms,
                         switchDistance=8.0*angstroms,
                         constraints=HBonds)

ref_atoms = [{'name': 'CA', 'resname': 'SER', 'resid': '192','chain':0},{'name': 'CA', 'resname': 'GLY', 'resid': '204','chain':0}]

# Define GCMC Sampler
gcmc_mover = grand.samplers.StandardGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=298*kelvin,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=5.5*angstroms,
                                                      log='1c5t-prod2.log',
                                                      excessChemicalPotential=-6.09*kilocalorie_per_mole,
                                                      standardVolume=30.345*angstroms**3,
                                                      ghostFile='1c5t-prod2.txt',
                                                      dcd='1c5t-prod2.dcd',
                                                      rst='1c5t-prod2.rst7',
                                                      overwrite=True)


# Define integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picosecond)

# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object
simulation = Simulation(pdb.topology, system, integrator, platform)

# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
#simulation.context.setVelocitiesToTemperature(298*kelvin)
simulation.context.setPositions(rst7.getPositions())
simulation.context.setVelocities(rst7.getVelocities())
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Prepare the GCMC sphere
gcmc_mover.initialise(simulation.context, ghosts[-1])


# Add StateDataReporter for production
simulation.reporters.append(StateDataReporter(stdout,
                                              1000,
                                              step=True,
                                              potentialEnergy=True,
                                              temperature=True,
                                              volume=True,
                                              speed=True))

# Run GCMC/MD equilibration (500k GCMC moves over 10 ns - 50 moves every 1 ps)
for i in range(2500):
    simulation.step(500)
    gcmc_mover.move(simulation.context, 50)
    # Write out a frame every 20 ps
#    if (i+1) % 20 == 0:
    gcmc_mover.report(simulation)


