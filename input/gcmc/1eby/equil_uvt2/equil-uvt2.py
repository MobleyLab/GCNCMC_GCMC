from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
from sys import stdout
import numpy as np
import mdtraj
import grand

# Load PDB
pdb = PDBFile('1eby-npt.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='1eby-ghosts2.pdb')

# Create system

ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')

system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=10.0*angstroms,
                         switchDistance=8.0*angstroms,
                         constraints=HBonds)

# Define reference atoms for the GCMC sphere
ref_atoms = [{'name': 'CA', 'resname': 'GLY', 'resid': '49','chain':0},{'name': 'CA', 'resname': 'GLY', 'resid': '49','chain':1}]

# Define GCMC Sampler
gcmc_mover = grand.samplers.StandardGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=277*kelvin,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=4.2*angstroms,
                                                      log='1eby-equil-uvt2.log',
                                                      excessChemicalPotential=-6.34*kilocalorie_per_mole,
                                                      standardVolume=29.823*angstroms**3,
                                                      ghostFile='1eby-uvt2-ghosts.txt',
                                                      dcd='1eby-equil-uvt2.dcd',
                                                      rst='1eby-equil-uvt2.rst7',
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

# Run GCMC/MD equilibration (100k GCMC moves over 500 ps - 200 moves every 1 ps)
for i in range(500):
#for i in range(10):
    print(i)
    simulation.step(500)
    gcmc_mover.move(simulation.context, 200)
    gcmc_mover.report(simulation)

# Remove ghosts and write out a PDB
ghost_resids = gcmc_mover.getWaterStatusResids(0)
positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
pdb.topology, pdb.positions = grand.utils.remove_ghosts(pdb.topology, positions,
                                                        ghosts=ghost_resids,
                                                        pdb='1eby-uvt2.pdb')

