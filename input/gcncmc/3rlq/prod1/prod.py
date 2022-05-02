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
pdb = PDBFile('3rlq-uvt2.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='3rlq-ghosts.pdb')


ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=10.0*angstroms,
                         switchDistance=8.0*angstroms,
                         constraints=HBonds)

ref_atoms = [{'name': 'CA', 'resname': 'SER', 'resid': '44', 'chain': 0},{'name': 'CA', 'resname': 'LYS', 'resid': '177', 'chain': 0}]

# Define integrator
integrator = BAOABIntegrator(286*kelvin, 1.0/picosecond, 0.002*picosecond)
# Define GCMC Sampler
# Define GCMC Sampler
gcncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=286*kelvin,
                                                      integrator=integrator,
                                                      nPertSteps=99, nPropStepsPerPert=40,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=6.5*angstroms,
                                                      log='3rlq-prod.log',
                                                      excessChemicalPotential=-6.19*kilocalorie_per_mole,
                                                      standardVolume=30.035*angstroms**3,
                                                      ghostFile='3rlq-prod.txt',
                                                      dcd='3rlq-prod.dcd',
                                                      rst='3rlq-prod.rst7',
                                                      overwrite=True)


# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object
simulation = Simulation(pdb.topology, system, gcncmc_mover.compound_integrator, platform)
# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(286*kelvin)
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

#simulation.reporters.append(DCDReporter('3rlq-prod.dcd', 100))  
# Run GCMC/MD equilibration (500k GCMC moves over 10 ns - 50 moves every 1 ps)
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

sys.exit()

#
# Format trajectory for visualisation
#

# Remove ghost waters from GCMC region
trj = grand.utils.shift_ghost_waters(ghost_file='3rlq-prod.txt',
                                     topology='3rlq-ghosts.pdb',
                                     trajectory='3rlq-prod.dcd')

# Centre the trajectory on a particular residue
trj = grand.utils.recentre_traj(t=trj, resname='ASN', name='CA', resid=36)

# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output='3rlq-gcmc-md.dcd')

# Write out a PDB trajectory of the GCMC sphere
grand.utils.write_sphere_traj(radius=5.5,
                              ref_atoms=ref_atoms,
                              topology='3rlq-ghosts.pdb',
                              trajectory='3rlq-gcmc-md.dcd',
                              output='gcmc_sphere.pdb',
                              initial_frame=True)

# Cluster water sites
grand.utils.cluster_waters(topology='3rlq-ghosts.pdb',
                           trajectory='3rlq-gcmc-md.dcd',
                           sphere_radius=5.5,
                           ref_atoms=ref_atoms,
                           cutoff=2.4,
                           output='watclusts.pdb')


