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
pdb = PDBFile('4zlz-uvt2.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='4zlz-ghosts.pdb')


ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=10.0*angstroms,
                         switchDistance=8.0*angstroms,
                         constraints=HBonds)

ref_atoms = [{'name': 'CA', 'resname': 'GLY', 'resid': '23', 'chain': 0},{'name': 'CA', 'resname': 'ILE', 'resid': '85', 'chain': 0}]

# Define integrator
integrator = BAOABIntegrator(277*kelvin, 1.0/picosecond, 0.002*picosecond)
# Define GCMC Sampler
gcncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=277*kelvin,
                                                      integrator=integrator,
                                                      nPertSteps=99, nPropStepsPerPert=40,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=5.5*angstroms,
                                                      log='4zlz-prod.log',
                                                      excessChemicalPotential=-6.34*kilocalorie_per_mole,
                                                      standardVolume=29.823*angstroms**3,
                                                      ghostFile='4zlz-prod.txt',
                                                      dcd='4zlz-prod.dcd',
                                                      rst='4zlz-prod.rst7',
                                                      overwrite=True)



# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')


# Create simulation object
simulation = Simulation(pdb.topology, system, gcncmc_mover.compound_integrator, platform)
# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(277*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Prepare the GCMC sphere
gcncmc_mover.initialise(simulation.context, ghosts)


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

#simulation.reporters.append(DCDReporter('2qbs-prod.dcd', 100))  
# Run GCMC/MD equilibration (500k GCMC moves over 10 ns - 50 moves every 1 ps)
for i in range(500):
    simulation.step(500)
    gcncmc_mover.move(simulation.context, 1)
#    gcmc_mover.move(simulation.context, 50)
    # Write out a frame every 20 ps
#    if (i+1) % 20 == 0:
    gcncmc_mover.report(simulation)


