import os, shutil, subprocess
import simtk.unit as unit
import simtk.openmm.app as app
import simtk.openmm as openmm

from paprika.restraints.utils import parse_window
from paprika.restraints.openmm import apply_positional_restraints, apply_dat_restraint
from paprika.build.system import TLeap
from paprika.build import align, dummy
from paprika import restraints, analysis, utils
from paprika.io import load_restraints, save_restraints

import numpy as np
import parmed as pmd

def APR_build(window_list, work_dir, guest_id, guest, solvate=False):
    base_name = "complex"
    complex_dir = "complex"
    for window in window_list:
        # Moving files in window's directory
        folder = os.path.join(work_dir, window)
        if not os.path.isdir(folder):
            os.makedirs(os.path.join(work_dir, window))
        if window[0] == "a":
            shutil.copy(os.path.join(complex_dir, f"{base_name}.prmtop"),
                        os.path.join(work_dir, window, f"{base_name}.prmtop"))
            shutil.copy(os.path.join(complex_dir, f"{base_name}.rst7"),
                        os.path.join(work_dir, window, f"{base_name}.rst7"))
            shutil.copy(os.path.join(complex_dir, f"{base_name}.pdb"),
                        os.path.join(work_dir, window, f"{base_name}.pdb"))

        elif window[0] == 'r':
            shutil.copy(os.path.join(work_dir, 'p039', f"{base_name}.prmtop"),
                        os.path.join(work_dir, window, f"{base_name}.prmtop"))
            shutil.copy(os.path.join(work_dir, 'p039', f"{base_name}.rst7"),
                        os.path.join(work_dir, window, f"{base_name}.rst7"))
            shutil.copy(os.path.join(work_dir, 'p039', f"{base_name}.pdb"),
                        os.path.join(work_dir, window, f"{base_name}.pdb"))
            
        elif window[0] == "p":
            structure = pmd.load_file(os.path.join(complex_dir, f"{base_name}.prmtop"),
                                      os.path.join(complex_dir, f"{base_name}.rst7"), structure = True)
            target_difference = guest_restraints[0].phase['pull']['targets'][int(window[1:])] -\
                                guest_restraints[0].pull['target_initial']

            for atom in structure.atoms:
                if atom.residue.name == guest_id:
                    atom.xz += target_difference

            structure.save(os.path.join(work_dir, window, f"{base_name}.prmtop"), overwrite=True)
            structure.save(os.path.join(work_dir, window, f"{base_name}.rst7"), overwrite=True)
            structure.write_pdb(os.path.join(work_dir, window, f"{base_name}.pdb"), renumber=False, write_links=False)
    
        # Current window
        window_number, phase = parse_window(window)
        print(f"Creating XML for in window {window}")

        #Solvate AMBER
        if solvate is True:
            system = TLeap()
            system.output_prefix = base_name
            system.output_path = folder
            system.target_waters = 2210
            system.neutralize = True
            system.add_ions = ['Na+', 6, 'Cl-', 6]
            system.pbc_type = 'rectangular'
            system.template_lines = ['source leaprc.GLYCAM_06j-1',
                                     'source leaprc.gaff2',
                                     #'source leaprc.protein.ff19SB',
                                     'source leaprc.water.tip3p',
                                     f'loadamberparams ../../{complex_dir}/dummy.frcmod',
                                     f'loadamberparams ../../{complex_dir}/{guest}.frcmod', 
                                     f'loadoff ../../{complex_dir}/{guest}.lib',
                                     f'loadoff ../../{complex_dir}/dm1.lib',
                                     f'loadoff ../../{complex_dir}/dm2.lib', 
                                     f'loadoff ../../{complex_dir}/dm3.lib', 
                                     f'model = loadpdb system.pdb']
            system.build(clean_files=False)
            system.repartition_hydrogen_mass()
        
        # Load Amber
        prmtop = app.AmberPrmtopFile(os.path.join(folder, f'{base_name}.prmtop'))
        inpcrd = app.AmberInpcrdFile(os.path.join(folder, f'{base_name}.rst7'))

        # Create PDB file
        with open(os.path.join(folder, 'system.pdb'), 'w') as file:
            app.PDBFile.writeFile(prmtop.topology, inpcrd.positions, file, keepIds=True)

        # Create an OpenMM system from the Amber topology
        if solvate is True:          
            system = prmtop.createSystem(nonbondedMethod=app.PME, constraints=app.HBonds,
                                         nonbondedCutoff=9 * unit.angstroms,
                                         hydrogenMass=3.024 * unit.daltons)

        else:
            system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
                                         implicitSolvent=app.HCT)

        # Apply positional restraints on the dummy atoms
        apply_positional_restraints(os.path.join(folder, 'system.pdb'), system, force_group=15)

        # Apply host static restraints
        for restraint in static_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=10)

        # Apply guest restraints
        for restraint in guest_restraints:
            apply_dat_restraint(system, restraint, phase, window_number, force_group=11)

        # Save OpenMM system to XML file
        system_xml = openmm.XmlSerializer.serialize(system)
        with open(os.path.join(folder, 'system.xml'), 'w') as file:
            file.write(system_xml)
