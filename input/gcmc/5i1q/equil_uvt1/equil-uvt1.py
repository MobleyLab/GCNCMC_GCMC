import sys
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
from sys import stdout
import numpy as np
import mdtraj
import grand

# Load PDB
pdb = PDBFile('min_complex.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='5i1q-ghosts.pdb')

ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=10.0*angstroms,
                         switchDistance=8.0*angstroms,
                         constraints=HBonds)

ref_atoms = [{'name': 'CA', 'resname': 'GLU', 'resid': '51', 'chain': 0},{'name': 'CA', 'resname': 'ASN', 'resid': '83', 'chain': 0}]

# Define GCMC Sampler
gcmc_mover = grand.samplers.StandardGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=277*kelvin,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=7.0*angstroms,
                                                      log='5i1q-equil-uvt1.log',
                                                      excessChemicalPotential=-6.34*kilocalorie_per_mole,
                                                      standardVolume=29.823*angstroms**3,
                                                      ghostFile='5i1q-uvt1-ghosts.txt',
                                                      dcd='5i1q-equil-uvt1.dcd',
                                                      rst='5i1q-equil-uvt1.rst7',
                                                      overwrite=True)

# Define integrator
integrator = BAOABIntegrator(277*kelvin, 1.0/picosecond, 0.002*picosecond)
# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object
simulation = Simulation(pdb.topology, system, integrator, platform)
# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(277*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Prepare the GCMC sphere
gcmc_mover.initialise(simulation.context, ghosts)
# Remove all waters currently in the sphere to reduce bias
gcmc_mover.deleteWatersInGCMCSphere()

# Start with 10k moves
#gcmc_mover.move(simulation.context, 10000)
print('Start with 10k moves')
gcmc_mover.move(simulation.context, 10000)


# Run GCMC/MD equilibration (100k GCMC moves over 1 ps - 1000 moves every 10 fs)
#for i in range(100):
for i in range(100):
    print(i)
    gcmc_mover.move(simulation.context, 1000)
    gcmc_mover.report(simulation)
    simulation.step(5)

# Remove ghosts and write out a PDB
# Remove ghosts and write out a PDB
ghost_resids = gcmc_mover.getWaterStatusResids(0)
positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
pdb.topology, pdb.positions = grand.utils.remove_ghosts(pdb.topology, positions,
                                                        ghosts=ghost_resids,
                                                        pdb='5i1q-uvt1.pdb')


