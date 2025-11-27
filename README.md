# Attach-Pull-Release (APR) Simulations of β-Cyclodextrin

The Attach-Pull-Release (APR) method is a molecular dynamics (MD) based free energy calculation technique used to estimate binding affinities in host–guest systems. It works by restraining, pulling, and releasing the guest molecule from the host, allowing accurate computation of binding thermodynamics:

    Attach: The guest molecule is gradually restrained to the host binding site using harmonic restraints. This ensures the system is well-defined and avoids sampling issues.

    Pull: The guest is then pulled away from the host along a defined reaction coordinate, typically using umbrella sampling windows.

    Release: Finally, restraints are removed so the guest is free in bulk solvent. This step accounts for the entropic contribution of binding.

Together, these steps yield the absolute binding free energy by integrating over the restraint work and the potential of mean force along the pulling coordinate.

> ### Note
> 
> This repository was built in 2021 and is no longer actively maintained
>


