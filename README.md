# Attach-Pull-Release (APR) Simulations of β-Cyclodextrin

## Background
The Attach-Pull-Release (APR) method is a molecular dynamics (MD) based free energy calculation technique used to estimate binding affinities in host–guest systems. It works by restraining, pulling, and releasing the guest molecule from the host, allowing accurate computation of binding thermodynamics:

__Attach__: The guest molecule is gradually restrained to the host binding site using harmonic restraints. This ensures the system is well-defined and avoids sampling issues.
__Pull__: The guest is then pulled away from the host along a defined reaction coordinate, typically using umbrella sampling windows.
__Release__: Finally, restraints are removed so the guest is free in bulk solvent. This step accounts for the entropic contribution of binding.

Together, these steps yield the absolute binding free energy by integrating over the restraint work and the potential of mean force along the pulling coordinate.

> ### Note
> 
> This repository was built in 2021 and is no longer actively maintained.
>


## Overview
The code in `/src/bcdmd` is structured more like a collection of scripts than a roboust python package, wrapping functionalities of paprika, ambertools and openmm into a simple workflow for bCD host-guest complexes.  As their file names imply, four key tasks are addressed:

| file | purpose|
|------|--------|
| _build.py_ | executes tleap to generate the force-field parameters for the host-guest complex while positioning dummy-atoms needed for the APR setup |
| _simbuilder.py_ | builds the simulations in openmm format for the guest molecules gradually pulled out of the host's binding pocket, while adding explicit solvent molecules and ions |
| _simulate.py_ | executes the openmm simulations of each window, initialising or appending to existing simulation data |
| _analysis.py_ | executes the free-energy calculation for the guest binding in the host based on the collected MD trajectories |