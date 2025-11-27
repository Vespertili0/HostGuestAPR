# Simple APR Simulation Scripts for β-Cyclodextrin Host–Guest Binding

## Background
The Attach-Pull-Release (APR) method is a molecular dynamics (MD) based free energy calculation technique used to estimate binding affinities in host–guest systems. It works by restraining, pulling, and releasing the guest molecule from the host, allowing accurate computation of binding thermodynamics:

1. __Attach__: The guest molecule is gradually restrained to the host binding site using harmonic restraints. This ensures the system is well-defined and avoids sampling issues.
2. __Pull__: The guest is then pulled away from the host along a defined reaction coordinate, typically using umbrella sampling windows.
3. __Release__: Finally, restraints are removed so the guest is free in bulk solvent. This step accounts for the entropic contribution of binding.

![](notebooks/files/apr_graphic.png)

Together, these steps yield the absolute binding free energy by integrating over the restraint work and the potential of mean force along the pulling coordinate.

## Overview

> ### :warning: Note
> 
> This repository was built in 2021 and is no longer actively maintained.
>

The code in `/src/bcdmd` is structured more like a collection of scripts than a roboust python package, wrapping functionalities of [paprika](https://github.com/GilsonLabUCSD/pAPRika), [ambertools](https://ambermd.org/AmberTools.php) and [openmm](https://github.com/openmm/openmm) into a simple workflow for β-CD host-guest complexes. As their file names imply, four key tasks are addressed:

| file | purpose|
|------|--------|
| _build.py_ | executes tleap to generate the force-field parameters for the host-guest complex while positioning dummy-atoms needed for the APR setup |
| _simbuilder.py_ | builds the simulations in openmm format for the guest molecules gradually pulled out of the host's binding pocket, while adding explicit solvent molecules and ions |
| _simulate.py_ | executes the openmm simulations of each window, initialising or appending to existing simulation data |
| _analysis.py_ | executes the free-energy calculation for the guest binding in the host based on the collected MD trajectories |

### Usage
Below is a typical workflow using the scripts in `src/bcdmd` for β-CD host–guest APR simulations.

Each step can be run independently, allowing flexible and modular execution.

#### 1. Build the Host–Guest System
Here, the General AMBER Force Field (GAFF) is used. Alternatively, GLYCAM_06j-1 is also implemented by using _Glycam_ as keyword instead. In this case, the β-CD in the _complex.pdb_ has to be written in the glycam-specific format.
First, the force-field parameters are applied, then the complex is aligned in z-direction of the periodic cell and the dummy-atoms are added to anchor the complex. Hydrogen Mass Repartitioning (HMR) is applied to enable larger time steps (4 fs) during MD simulations.

```python
from bcdmd.build import HostGuestComplexSetup

apr_dir = "/path/to/APR_sims/"
guest_id = "ANA" #ligand id used in complex.pdb file
complex_pdb = "complex.pdb"
hgcs = HostGuestComplexSetup(apr_dir, complex_pdb, guest_id, "GAFF2")

hgcs.parameterise_structure()
hgcs.positioning_complex(":ANA@C9", ":ANA@C3")
hgcs.provide_dummies()
```

#### 2. Prepare Simulation Windows
The simulation windows are build, gradually pulling the guest molecule along the z-axis (here C9 -> C3) out of the guest molecule. Each window is written to a separate directory with its corresponding structure and force-field files for execution with openmm. The existing windows are solvated during a second build-process.

```python
from bcdmd.simbuilder import HostGuestComplexSimulationBuilder

hgcsb = HostGuestComplexSimulationBuilder(apr_dir, True, "GAFF")
hgcsb.apply_restraints(":ANA@C9", ":ANA@C3")
# first complex only
hgcsb.APR_build(guest_id=pdb_id, guest=pdb_id, solvate=False)
# then add solvation
hgcsb.APR_build(guest_id=pdb_id, guest=pdb_id, solvate=True)
```

#### 3. Run Simulations
Simulations are executed via openmm, applying pre-set conditions of 300K with a Langevin-Middle integrator. Initial runs execute a minimisation and NPT sequence automatically. 

```python
from bcdmd.simulate import MDrunAPR

# Run MD for each window
# use restart=False for initial runs
# and restart=True to continue
FRAME_list = ["a000"]
loop = 0
while FRAME_list and loop < 5:
    loop += 1
    FRAME_list = MDrunAPR(FRAME_list, f"{apr_dir}/apr_windows/")
    print(FRAME_list, f"...after: {loop} loop")
```

#### 4. Analyze Results
Wraps the analysis steps as implemented in [_paprika_](https://paprika.readthedocs.io/en/latest/workflow.html#analysis).

```python
from bcdmd.analysis import APRanalysis

fe = APRanalysis(apr_dir, pdb_id)
fe.pull_Data()
fe.run_Thermo()
fe.summarize_APR()
fe.plot_SEMatrix("pull")
```

> See [`notebooks/APR_system_GAFF.ipynb`](notebooks/APR_system_GAFF.ipynb) for a full workflow example.

