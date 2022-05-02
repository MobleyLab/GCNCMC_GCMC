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
pdb = PDBFile('5i1q_ghosts_gcmc.pdb')

# Load restart file
rst7 = AmberInpcrdFile('5i1q_prod5_gcmc.rst7')

# Shouldn't need to add ghosts as these can just be read in from before (all frames contained below)
ghosts = grand.utils.read_ghosts_from_file('5i1q_prod5_gcmc.txt')

ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=10.0*angstroms,
                         switchDistance=8.0*angstroms,
                         constraints=HBonds)

ref_atoms = [{'name': 'CA', 'resname': 'GLU', 'resid': '51', 'chain': 0},{'name': 'CA', 'resname': 'ASN', 'resid': '83', 'chain': 0}]

# Define integrator
integrator = BAOABIntegrator(277*kelvin, 1.0/picosecond, 0.002*picosecond)


# Define GCMC Sampler

gcncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=277*kelvin,
                                                      integrator=integrator,
                                                      nPertSteps=99, nPropStepsPerPert=40,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=7.0*angstroms,
                                                      log='5i1q-prod1.log',
                                                      excessChemicalPotential=-6.19*kilocalorie_per_mole,
                                                      standardVolume=29.823*angstroms**3,
                                                      ghostFile='5i1q-prod1.txt',
                                                      dcd='5i1q-prod1.dcd',
                                                      rst='5i1q-prod1.rst7',
                                                      overwrite=True)



# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object
simulation = Simulation(pdb.topology, system, gcncmc_mover.compound_integrator, platform)

# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
#simulation.context.setVelocitiesToTemperature(278*kelvin)
simulation.context.setPositions(rst7.getPositions())
simulation.context.setVelocities(rst7.getVelocities())
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Prepare the GCMC sphere
gcncmc_mover.initialise(simulation.context, ghosts[-1])


# Add StateDataReporter for production
simulation.reporters.append(StateDataReporter(stdout,
                                              1000,
                                              step=True,
                                              potentialEnergy=True,
                                              temperature=True,
                                              volume=True,
                                              speed=True))

# Reset GCMC statistics
gcncmc_mover.reset()


# Run GCMC/MD equilibration (500k GCMC moves over 10 ns - 50 moves every 1 ps)
for i in range(500):
    simulation.step(500)
    gcncmc_mover.move(simulation.context, simulation, 1)
#    gcmc_mover.move(simulation.context, 50)
    # Write out a frame every 20 ps
#    if (i+1) % 20 == 0:
    gcncmc_mover.report(simulation)



