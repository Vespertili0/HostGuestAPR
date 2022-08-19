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

class HostGuestComplexSetup:
    def __init__(self, work_dir, cx_name, guest_id, FF):
        self.work_dir = work_dir
        self.cx_name = cx_name
        self.guest_id = guest_id
        self.FF = FF
    
    def _run_tleap_prep(self, name, template_list):
        system = TLeap()
        system.output_path = 'complex/'
        system.output_prefix = name
        system.pbc_type = None
        system.neutralize = False
        system.template_lines = template_list
        system.build(clean_files=False)
        system.run()
        print('DONE...check tleap-log')

    def parameterise_structure(self):
        '''
        Parameters
        ----------
        cx_name : str
            name of pdb-file w/o .pdb-ending
        guest : str
            pdb-ID of guest-molecule in complex
        FF : str
            Glycam or GAFF2
        '''
        if self.FF == 'Glycam':
            FF_string = 'source leaprc.GLYCAM_06j-1'
        if self.FF == 'GAFF2': 
            FF_string = 'bCD = loadmol2 bCD.mol2\nloadamberparams bCD.frcmod'#\nloadoff bCD.lib'
        folder = 'Docked_Structures/'
        template_list = ['source leaprc.gaff2',
                         FF_string,
                         #'source leaprc.protein.ff19SB',
                         f'{self.guest_id} = loadmol2 {self.guest_id}.mol2',
                         f'loadamberparams {self.guest_id}.frcmod',
                       #  f'loadoff {self.guest_id}.lib',
                         f'complex = loadpdb {self.cx_name}.pdb',
                        # f'saveamberparm complex {cx_name}.prmtop {cx_name}.rst7'
                        ]
        self._run_tleap_prep('vac', template_list)

    def positioning_complex(self, G1, G2):
        folder = os.path.join(self.work_dir, 'complex')
        structure = pmd.load_file(f'{folder}/vac.prmtop',
                                  f'{folder}/vac.rst7', structure=True)
        structure = align.zalign(structure, G1, G2)
        structure = dummy.add_dummy(structure, residue_name='DM1', z=-6)
        structure = dummy.add_dummy(structure, residue_name='DM2', z=-9)
        structure = dummy.add_dummy(structure, residue_name='DM3', z=-11.2, y=2.2)
        structure.save(f'{folder}/aligned_dum.prmtop', overwrite=True)
        structure.save(f'{folder}/aligned_dum.rst7', overwrite=True)
        #structure.write_pdb(f'{folder}/aligned_dum.pdb', renumber=False, write_links=True)
        structure.save(f'{folder}/aligned_dum.pdb', overwrite=True)
        dummy.write_dummy_frcmod(filepath=f'{folder}/dummy.frcmod')
        for i in range(1,4):
            dm_name = folder + f'/dm{i}.mol2'
            dummy.write_dummy_mol2(residue_name=f'DM{i}', filepath=dm_name)
            
    def _modify_aligned_pdb(self, file_loc):
        '''
        Relevant in any FF
        [x] rename dummy-atoms from Pb to Du
        
        Relevant for usage of Glycam as FF
        [] convert with openbabel to add CONECT section
        [] add LINK-section for 4GA-unit linker
        [] check CONECT-section for closed bCD-ring
        [] remove DM1-DM2-DM3 from pdb-CONECT-section (last three entries)
        '''
        pdb = open(f'{file_loc}/aligned_dum.pdb').read()
        pdb = pdb.replace('  PB', '  Du')
        with open(f'{file_loc}/aligned_dum.pdb', 'w') as file:
            file.write(pdb)

    def provide_dummies(self):
        folder = os.path.join(self.work_dir, 'complex')
        self._modify_aligned_pdb(folder)
        if self.FF == 'GAFF2':
            template_list = ['source leaprc.gaff2',
                            'bCD = loadmol2 bCD.mol2',
                            'loadamberparams bCD.frcmod']
        if self.FF == 'Glycam':
            template_list = ['source leaprc.GLYCAM_06j-1','source leaprc.gaff2']
        template_rest = [f'loadamberparams dummy.frcmod',
                         f'{self.guest_id} = loadmol2 {self.guest_id}.mol2',
                         f'loadamberparams {self.guest_id}.frcmod',
                        # f'loadoff {self.guest_id}.lib',
                         f'DM1 = loadmol2 dm1.mol2',
                         f'DM2 = loadmol2 dm2.mol2',
                         f'DM3 = loadmol2 dm3.mol2',
                         f'saveoff DM1 dm1.lib',
                         f'saveoff DM2 dm2.lib',
                         f'saveoff DM3 dm3.lib',
                         f'complex = loadpdb aligned_dum.pdb',
                         'check complex']
        template_list.extend(template_rest)
        self._run_tleap_prep('complex', template_list)
        shutil.copy(os.path.join(folder, 'aligned_dum.pdb'), os.path.join(folder, 'complex.pdb'))