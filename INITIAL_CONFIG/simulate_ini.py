#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:46:30 2023

@author: albert
"""
import simtk.openmm as openmm
from simtk.openmm import unit
from simtk.openmm import app
from simtk.openmm import XmlSerializer
from analyse_ini import *
import time
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--seq',nargs='?',const='', type=str)
parser.add_argument('--temp',nargs='?',const='', type=int)
parser.add_argument('--cutoff',nargs='?',const='', type=float)
parser.add_argument('--windows',nargs='?',const='', type=int)
parser.add_argument('--n_chains',nargs='?',const='', type=int)
parser.add_argument('--steps',nargs='?',const='', type=int)
parser.add_argument('--L',nargs='?',const='', type=float)
args = parser.parse_args()

def simulate(residues,name,prot,temp,cutoff,pdb_file,steps,pdb_file2,n_clust,perc,n_chains, L):
    residues = residues.set_index('one')

    lj_eps, fasta, types, MWs = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa = genParamsDH(residues,name,prot,temp)

    N = len(fasta)
    # set parameters

    margin=8.
    marginz=0.38*N/2+8.
    Nsteps=int(6e7)

    system = openmm.System()

    # set box vectors
    a = unit.Quantity(np.zeros([3]), unit.nanometers)
    a[0] = L * unit.nanometers
    b = unit.Quantity(np.zeros([3]), unit.nanometers)
    b[1] = L * unit.nanometers
    c = unit.Quantity(np.zeros([3]), unit.nanometers)
    c[2] = L * unit.nanometers
    system.setDefaultPeriodicBoxVectors(a, b, c)

    pdb = app.pdbfile.PDBFile(pdb_file)

    for chain in range(n_chains):
        system.addParticle((residues.loc[prot.fasta[0]].MW+2)*unit.amu)
        count=0
        for a in prot.fasta[1:-1]:
            system.addParticle(residues.loc[a].MW*unit.amu)
            if count==int(len(prot.fasta)/2) and chain > n_clust:
                system.setParticleMass(count, 0*unit.amu)
            count+=1

        system.addParticle((residues.loc[prot.fasta[-1]].MW+16)*unit.amu)

    hb = openmm.openmm.HarmonicBondForce()

    energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6-shift),4*eps*((s/r)^12-(s/r)^6-l*shift)+eps*(1-l))'
    ah = openmm.openmm.CustomNonbondedForce(energy_expression+'; s=0.5*(s1+s2); l=0.5*(l1+l2); shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6')

    ah.addGlobalParameter('eps',lj_eps*unit.kilojoules_per_mole)
    ah.addGlobalParameter('rc',cutoff*unit.nanometer)
    ah.addPerParticleParameter('s')
    ah.addPerParticleParameter('l')
    
    print('cutoff distance: ',cutoff*unit.nanometer)
 
    yu = openmm.openmm.CustomNonbondedForce('q*(exp(-kappa*r)/r-shift); q=q1*q2')
    yu.addGlobalParameter('kappa',yukawa_kappa/unit.nanometer)
    yu.addGlobalParameter('shift',np.exp(-yukawa_kappa*4.0)/4.0/unit.nanometer)
    yu.addPerParticleParameter('q')

    for j in range(n_chains):
        begin = j*N
        end = j*N+N
       
        for a,e in zip(prot.fasta,yukawa_eps):
            yu.addParticle([e*unit.nanometer*unit.kilojoules_per_mole])
            ah.addParticle([residues.loc[a].sigmas*unit.nanometer, residues.loc[a].lambdas*unit.dimensionless])

        for i in range(begin,end-1):
            hb.addBond(i, i+1, 0.38*unit.nanometer, 8033.0*unit.kilojoules_per_mole/(unit.nanometer**2))
            yu.addExclusion(i, i+1)
            ah.addExclusion(i, i+1)

    yu.setForceGroup(0)
    ah.setForceGroup(1)
    yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    hb.setUsesPeriodicBoundaryConditions(True)
    yu.setCutoffDistance(4*unit.nanometer)
    ah.setCutoffDistance(cutoff*unit.nanometer)

    system.addForce(hb)
    system.addForce(yu)
    system.addForce(ah)

    print('Using pbc in LJ: ',ah.usesPeriodicBoundaryConditions())
    print('Using pbc in Yukawa: ', yu.usesPeriodicBoundaryConditions())
    print('Using pbc in Harmonic bond force: ',hb.usesPeriodicBoundaryConditions())

    integrator = openmm.openmm.LangevinIntegrator(temp*unit.kelvin,0.01/unit.picosecond,0.01*unit.picosecond)

    print('Friction: ',integrator.getFriction(),'Temperature:',integrator.getTemperature())

    platform = openmm.Platform.getPlatformByName('CUDA')

    simulation = app.simulation.Simulation(pdb.topology, system, integrator, platform, dict(CudaPrecision='mixed'))

    check_point = name+'/{:d}/restart.chk'.format(temp)

    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    simulation.reporters.append(app.dcdreporter.DCDReporter(name+'/{:d}/{:s}.dcd'.format(temp,name),steps//10))
    # simulation.reporters.append(app.pdbreporter.PDBReporter('config/top_eq_0.50.pdb',10000,enforcePeriodicBox=True))
    
    print("~~~ STARTING SIMULATION ~~~")
    simulation.reporters.append(app.statedatareporter.StateDataReporter('{:s}_{:d}.log'.format(name,temp),int(steps/10),
             step=True,speed=True,elapsedTime=True,separator='\t'))

    # simulation.runForClockTime(1/60.*unit.hour, checkpointFile=check_point, checkpointInterval=2*unit.hour)
    simulation.step(steps)
    simulation.saveCheckpoint(check_point)
    
    print("~~~ STARTING ANALYSIS ~~~")

    genDCD(residues,name,prot,temp,n_chains,pdb_file2,perc)

residues = pd.read_csv('residues.csv').set_index('three',drop=False)
proteins = pd.read_pickle('proteins.pkl')
#steps = 4000000
print('Protein name: ',args.seq,'\n Temperature: ',args.temp)
print('Steps: ',args.steps)

t0 = time.time()
for i in range(args.windows):
    perc=i/(args.windows-1)
    n_clust = int(args.n_chains*perc)
    pdb_file = 'config/top_{:4.2f}.pdb'.format(perc)
    pdb_file2 = 'config/top_eq_{:4.2f}.pdb'.format(perc)
    print('Generating '+pdb_file2)
    simulate(residues,args.seq,proteins.loc[args.seq],args.temp,args.cutoff,
             pdb_file,args.steps,pdb_file2,n_clust,perc,args.n_chains,args.L)
print('Timing {:.3f}'.format(time.time()-t0))
