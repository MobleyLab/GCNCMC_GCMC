import sys, os
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
from sys import stdout
import numpy as np
import mdtraj
import grand
import parmed as pmd


# Load PDB
pdb = PDBFile('1ec1-uvt2.pdb')
structure = pmd.load_file('1ec1-uvt2.pdb')
selection = '(!@H=&!(:HOH|:CL,NA))'
mask = pmd.amber.AmberMask(structure, selection)
mask_idx = [i for i in mask.Selected()]


# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='1ec1-ghosts.pdb')


ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml', 'LIG_test.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=10.0*angstroms,
                         switchDistance=8.0*angstroms,
                         constraints=HBonds)

ref_atoms = [{'name': 'CA', 'resname': 'GLY', 'resid': '49','chain':0},{'name': 'CA', 'resname': 'GLY', 'resid': '49','chain':1}]

force = CustomExternalForce('k_restr*periodicdistance(x, y, z, x0, y0, z0)^2')
force.addGlobalParameter("k_restr", 10.0*kilocalories_per_mole/angstroms**2)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")
for i, atom_crd in enumerate(structure.positions):
    if i in mask_idx:
        force.addParticle(i, atom_crd.value_in_unit(nanometers))
system.addForce(force)


# Define GCMC Sampler
gcmc_mover = grand.samplers.StandardGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=277*kelvin,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=5.5*angstroms,
                                                      log='1ec1-prod.log',
                                                      excessChemicalPotential=-6.34*kilocalorie_per_mole,
                                                      standardVolume=29.823*angstroms**3,
                                                      ghostFile='1ec1-prod.txt',
                                                      dcd='1ec1-prod.dcd',
                                                      rst='1ec1-prod.rst7',
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


# Add StateDataReporter for production
simulation.reporters.append(StateDataReporter(stdout,
                                              1000,
                                              step=True,
                                              potentialEnergy=True,
                                              temperature=True,
                                              volume=True,
                                              speed=True))

#simulation.reporters.append(DCDReporter('1ec1-prod.dcd', 100))  
# Run GCMC/MD equilibration (500k GCMC moves over 10 ns - 50 moves every 1 ps)
for i in range(2500):
    simulation.step(500)
    gcmc_mover.move(simulation.context, 50)
    # Write out a frame every 20 ps
#    if (i+1) % 20 == 0:
    gcmc_mover.report(simulation)

sys.exit()

#
# Format trajectory for visualisation
#

# Remove ghost waters from GCMC region
trj = grand.utils.shift_ghost_waters(ghost_file='1ec1-prod.txt',
                                     topology='1ec1-ghosts.pdb',
                                     trajectory='1ec1-prod.dcd')

# Centre the trajectory on a particular residue
trj = grand.utils.recentre_traj(t=trj, resname='ASN', name='CA', resid=36)

# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output='1ec1-gcmc-md.dcd')

# Write out a PDB trajectory of the GCMC sphere
grand.utils.write_sphere_traj(radius=5.5,
                              ref_atoms=ref_atoms,
                              topology='1ec1-ghosts.pdb',
                              trajectory='1ec1-gcmc-md.dcd',
                              output='gcmc_sphere.pdb',
                              initial_frame=True)

# Cluster water sites
grand.utils.cluster_waters(topology='1ec1-ghosts.pdb',
                           trajectory='1ec1-gcmc-md.dcd',
                           sphere_radius=5.5,
                           ref_atoms=ref_atoms,
                           cutoff=2.4,
                           output='watclusts.pdb')


