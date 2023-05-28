import os, shutil, subprocess
import simtk.unit as unit
import simtk.openmm.app as app
import simtk.openmm as openmm

from paprika.restraints.utils import parse_window
from paprika.restraints.openmm import apply_positional_restraints, apply_dat_restraint
from paprika.build.system import TLeap
from paprika.build.system.utils import PBCBox
from paprika.build import align, dummy
from paprika import restraints, analysis, utils
from paprika.io import load_restraints, save_restraints

import MDAnalysis as mda
import numpy as np
import parmed as pmd

class HostGuestComplexSimulationBuilder:
    def __init__(self, work_dir, include_release=False, FF='Glycam'):
        """
        include_release  : (bool)
            if True, release-phase will be generated
        """
        self.FF = FF
        self.work_dir= work_dir
        self.structure = pmd.load_file('complex/complex.prmtop',
                                       'complex/complex.rst7', structure=True)
        sp = [0, 0.4, 0.8, 1.5, 2.5, 4, 6, 8.5, 12, 18, 25, 35, 50, 75, 100]
        self.attach_fractions = [i/100 for i in sp]
        self.release_fractions = list()
        initial_distance_to_dummy = 6
        max_pull_distance = 16
        pull_interval = 0.4
        self.pull_distances = np.arange(0 + initial_distance_to_dummy,
                                       max_pull_distance + initial_distance_to_dummy,
                                       pull_interval)
        if include_release is True:
            sp.reverse()
            self.release_fractions = [i/100 for i in sp]         
        self.windows = [len(self.attach_fractions),
                        len(self.pull_distances),
                        len(self.release_fractions)]
        print(f'initializing attach-pull-release with: {self.windows} windows')
        
    def _create_StaticHostRestraint(self, mask_list, f_const):
        r = restraints.static_DAT_restraint(restraint_mask_list=mask_list,
                                            force_constant=f_const,
                                            ref_structure=self.structure,
                                            num_window_list=self.windows, amber_index=False)
        return r

    def _create_GuestRestraint(self, mask_list, rest_type):
        """rest_type: 'translation' - requires mask_list of 2-dim
                      'rotation' - requires mask_list of 3-dim"""
        r = restraints.DAT_restraint()
        r.topology = self.structure
        r.auto_apr = True
        r.continuous_apr = True
        r.amber_index = False
        if rest_type == 'translation':
            r.mask1 = mask_list[0]
            r.mask2 = mask_list[1]
            r.attach["target"] = self.pull_distances[0]          # Angstroms
            r.attach["fraction_list"] = self.attach_fractions
            r.attach["fc_final"] = 5.0                      # kcal/mol/Angstroms**2
            r.pull["target_final"] = 22                   # Angstroms
            r.pull["num_windows"] = self.windows[1]
            r.release['target'] = 22
            r.release['fraction_list'] = [1.0] * len(self.attach_fractions)
            r.release['fc_final'] = 5.0

        if rest_type == 'rotation':
            r.mask1 = mask_list[0]
            r.mask2 = mask_list[1]
            r.mask3 = mask_list[2]
            r.attach["target"] = 180.0                      # Degrees
            r.attach["fraction_list"] = self.attach_fractions
            r.attach["fc_final"] = 100.0                    # kcal/mol/radian**2
            r.pull["target_final"] = 180.0                  # Degrees
            r.pull["num_windows"] = self.windows[1]
            r.release['target'] = 180.0
            r.release['fraction_list'] = [1.0] * len(self.attach_fractions)
            r.release['fc_final'] = 100
        r.initialize()
        return r

    def _create_HostConformationalRestraint(self, mask_list):
        """ bCD-ring-integrity  - requires mask_list of 4-dim """
        r = restraints.DAT_restraint()
        r.topology = self.structure
        r.auto_apr = True
        r.continuous_apr = True
        r.amber_index = False
        r.mask1 = mask_list[0]
        r.mask2 = mask_list[1]
        r.mask3 = mask_list[2]
        r.mask4 = mask_list[3]
        f = lambda x: 108 if 'O5' in x else -114
        r.attach["target"] = f(mask_list[0])            # Degrees
        r.attach["fraction_list"] = self.attach_fractions
        r.attach["fc_final"] = 100.0                    # kcal/mol/radian**2
        r.pull["target_final"] = f(mask_list[0])        # Degrees
        r.pull["num_windows"] = self.windows[1]
        r.release['target'] = f(mask_list[0])   
        r.release['fraction_list'] = self.release_fractions
        r.release['fc_final'] = 100.0
        r.initialize()
        return r
    
    def _build_host_conf_mask(self):
        conf_Mask = list()
        for i in range(1, 8):
            f = lambda x: 1 if x==8 else x
            conf_Mask.append([f":{f(i+1)},@O5", f":{f(i+1)},@C1", f":{i},@O4", f":{i},@C4"])
            conf_Mask.append([f":{f(i+1)},@C1", f":{i},@O4", f":{i},@C4", f":{i},@C5"])
        return conf_Mask    
    
    def _convert_Glyc2Gaff(self):
        mol_glyc = mda.Universe(f'../mol2-files/bCD_glyc06_centered.pdb')
        mol_gaff = mda.Universe(f'../mol2-files/bCD_gaff.pdb')
        glyc_atom_names = list()
        gaff_atom_names = list()
        for a_chain, a_name in zip(list(mol_glyc.atoms.resids), list(mol_glyc.atoms.names)):
            glyc_atom_names.append(f':{a_chain},@{a_name}')
        for a_name in list(mol_gaff.atoms.names):
            gaff_atom_names.append(f':bCD@{a_name}')
        self.DICT_glyc2gaff = dict(zip(glyc_atom_names, gaff_atom_names))        
    
    def apply_restraints(self, G1, G2):
        """
        creates window-list with static restraints (e.g. host orientation)
        and dynamic restraints (e.g. rotation/transtation of guest, conformation of host)
        
        While static restraints do not change during the whole APR process and 
        do not affect the free energy, dynamic restraints do.
        
        returns window-list        
        """
        H1, H2, H3 = ":1,@C1", ":3,@C1", ":5,@C1"
        D1, D2, D3 = ":DM1", ":DM2", ":DM3"
        conf_MASK = self._build_host_conf_mask()
        if self.FF == 'GAFF':
            print('converting atom-names to GAFF')
            self._convert_Glyc2Gaff()
            H1 = self.DICT_glyc2gaff.get(H1)
            H2 = self.DICT_glyc2gaff.get(H2)
            H3 = self.DICT_glyc2gaff.get(H3)
            modified_conf_MASK = list()
            for mask_list in modified_conf_MASK:
                modified_mask = [self.DICT_glyc2gaff.get(glyc) for glyc in mask_list]
                modified_conf_MASK.append(modified_mask)
            conf_MASK = modified_conf_MASK.copy()
    
    ### static restraints ###
        self.static_restraints = list()
        static_MASK = [[D1, H1], [D2, D1, H1], [D3, D2, D1, H1],
                       [D1, H1, H2], [D2, D1, H1, H2], [D1, H1, H2, H3]]
        static_FORCE = [5.0, 100.0, 100.0, 100.0, 100.0, 100.0]
        for mask, force in zip(static_MASK, static_FORCE):
            self.static_restraints.append(self._create_StaticHostRestraint(mask, force))
    
    ### dynamic restraints ###
        self.dynamic_restraints = list()
        trans_rot_MASK = [[D1, G1], [D2, D1, G1], [D1, G1, G2]]
        TRANS_ROT = ['translation', 'rotation', 'rotation']
        for mask, trans_rot in zip(trans_rot_MASK, TRANS_ROT):
            self.dynamic_restraints.append(self._create_GuestRestraint(mask, trans_rot))
        
        for mask in conf_MASK:
            self.dynamic_restraints.append(self._create_HostConformationalRestraint(mask))
        self.window_list = restraints.create_window_list(self.dynamic_restraints)    
        save_restraints((self.static_restraints + self.dynamic_restraints),
                        filepath=f"{self.work_dir}/restraints.json")
        print('restraints saved to json-file')

    def APR_build(self, guest_id, guest, solvate=False, apr_dir='apr_windows'):
        base_name = "complex"
        complex_dir = "complex"
        for window in self.window_list:
            # Moving files in window's directory
            folder = os.path.join(apr_dir, window)
            if not os.path.isdir(folder):
                os.makedirs(os.path.join(apr_dir, window))
            if window[0] == "a":
                shutil.copy(os.path.join(complex_dir, f"{base_name}.prmtop"),
                            os.path.join(apr_dir, window, f"{base_name}.prmtop"))
                shutil.copy(os.path.join(complex_dir, f"{base_name}.rst7"),
                            os.path.join(apr_dir, window, f"{base_name}.rst7"))
                shutil.copy(os.path.join(complex_dir, f"{base_name}.pdb"),
                            os.path.join(apr_dir, window, f"{base_name}.pdb"))

            elif window[0] == 'r':
                shutil.copy(os.path.join(apr_dir, 'p039', f"{base_name}.prmtop"),
                            os.path.join(apr_dir, window, f"{base_name}.prmtop"))
                shutil.copy(os.path.join(apr_dir, 'p039', f"{base_name}.rst7"),
                            os.path.join(apr_dir, window, f"{base_name}.rst7"))
                shutil.copy(os.path.join(apr_dir, 'p039', f"{base_name}.pdb"),
                            os.path.join(apr_dir, window, f"{base_name}.pdb"))

            elif window[0] == "p":
                self.structure = pmd.load_file(os.path.join(complex_dir, f"{base_name}.prmtop"),
                                               os.path.join(complex_dir, f"{base_name}.rst7"),
                                               structure = True)
                target_difference = self.dynamic_restraints[0].phase['pull']['targets'][int(window[1:])] -\
                                    self.dynamic_restraints[0].pull['target_initial']

                for atom in self.structure.atoms:
                    if atom.residue.name == guest_id:
                        atom.xz += target_difference

                self.structure.save(os.path.join(apr_dir, window, f"{base_name}.prmtop"), overwrite=True)
                self.structure.save(os.path.join(apr_dir, window, f"{base_name}.rst7"), overwrite=True)
                self.structure.write_pdb(os.path.join(apr_dir, window, f"{base_name}.pdb"), renumber=False, write_links=False)

            # Current window
            window_number, phase = parse_window(window)
            print(f"Creating XML for in window {window}")

            #Solvate AMBER
            if solvate is True:
                if self.FF == 'GAFF':
                    template_list = ['source leaprc.gaff2',
                                    f'bCD = loadmol2 ../../{complex_dir}/bCD.mol2',
                                    f'loadamberparams ../../{complex_dir}/bCD.frcmod']
                if self.FF == 'Glycam':
                    template_list = ['source leaprc.GLYCAM_06j-1','source leaprc.gaff2']
                template_rest = ['source leaprc.water.tip3p',
                                 f'{guest} = loadmol2 ../../{complex_dir}/{guest}.mol2',
                                 f'DM1 = loadmol2 ../../{complex_dir}/dm1.mol2',
                                 f'DM2 = loadmol2 ../../{complex_dir}/dm2.mol2',
                                 f'DM3 = loadmol2 ../../{complex_dir}/dm3.mol2',
                                 f'loadamberparams ../../{complex_dir}/dummy.frcmod',
                                 f'loadamberparams ../../{complex_dir}/{guest}.frcmod', 
                                # f'loadoff ../../{complex_dir}/{guest}.lib',
                                 f'loadoff ../../{complex_dir}/dm1.lib',
                                 f'loadoff ../../{complex_dir}/dm2.lib', 
                                 f'loadoff ../../{complex_dir}/dm3.lib', 
                                 f'model = loadpdb system.pdb']
                template_list.extend(template_rest)
                
                system = TLeap()
                system.output_prefix = base_name
                system.output_path = folder
                system.target_waters = 2210
                system.neutralize = True
                system.add_ions = ['Na+', 6, 'Cl-', 6]
                system.pbc_type = PBCBox.rectangular
                system.template_lines = template_list
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
            for restraint in self.static_restraints:
                apply_dat_restraint(system, restraint, phase, window_number, force_group=10)

            # Apply guest restraints
            for restraint in self.dynamic_restraints:
                apply_dat_restraint(system, restraint, phase, window_number, force_group=11)

            # Save OpenMM system to XML file
            system_xml = openmm.XmlSerializer.serialize(system)
            with open(os.path.join(folder, 'system.xml'), 'w') as file:
                file.write(system_xml)
