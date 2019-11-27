#!/bin/bash

## Implosion with B fields
	
#SBATCH --job-name=implode_test

## Allocation Definition

#SBATCH --account=astro
#SBATCH --partition=astro

## Resources

## Nodes

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28

## Walltime - shouldn't be substantially longer than job needs

#SBATCH --time=1:00:00

## Memory per node - important to specify since default is small

## mox memory may be >=120G depending on node choice

#SBATCH --mem=100G

## Specify working directory for job

#SBATCH --chdir=/usr/lusers/sharpc42

## Email notification

#SBATCH --mail-type=ALL

#SBATCH --mail-user=sharpc42@uw.edu

#SBATCH --error=vtj-%j.err
#SBATCH --output=vtj-%j.out

## export all environmental variables to batch job session

#SBATCH --export=all

cd /gscratch/astro/sharpc42/athena-public-version/bin
module load icc_19-ompi_4.0.1
module load icc_19
mpirun athena -i athinput.blast -t 0:58:00 -d /gscratch/astro/sharpc42/Research/athena_runs/Hydrostatic/20191127_box-SH3_bcOOR_rad10_dt50_mb32_2D_cool
