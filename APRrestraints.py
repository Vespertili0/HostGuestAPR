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

def static_host_restraint(mask_list, f_const, structure, windows):
    r = restraints.static_DAT_restraint(restraint_mask_list=mask_list, force_constant=f_const,
                                        ref_structure=structure,
                                        num_window_list=windows, amber_index=False)
    return r

def guest_restraint(mask_list, rest_type, pull_distances, attach_fractions, windows):
    """rest_type: 'translation' - requires mask_list of 2-dim
                  'rotation' - requires mask_list of 3-dim"""
    r = restraints.DAT_restraint()
    r.topology = structure
    r.auto_apr = True
    r.continuous_apr = True
    r.amber_index = False
    if rest_type == 'translation':
        r.mask1 = mask_list[0]
        r.mask2 = mask_list[1]
        r.attach["target"] = pull_distances[0]          # Angstroms
        r.attach["fraction_list"] = attach_fractions
        r.attach["fc_final"] = 5.0                      # kcal/mol/Angstroms**2
        r.pull["target_final"] = 22                   # Angstroms
        r.pull["num_windows"] = windows[1]
        
        r.release['target'] = 22
        r.release['fraction_list'] = [1.0] * len(attach_fractions)
        r.release['fc_final'] = 5.0
        
    if rest_type == 'rotation':
        r.mask1 = mask_list[0]
        r.mask2 = mask_list[1]
        r.mask3 = mask_list[2]
        r.attach["target"] = 180.0                      # Degrees
        r.attach["fraction_list"] = attach_fractions
        r.attach["fc_final"] = 100.0                    # kcal/mol/radian**2
        r.pull["target_final"] = 180.0                  # Degrees
        r.pull["num_windows"] = windows[1]
        
        r.release['target'] = 180.0
        r.release['fraction_list'] = [1.0] * len(attach_fractions)
        r.release['fc_final'] = 100

    r.initialize()
    return r
  
  def conformational_host_restraint(mask_list, attach_fractions, release_fractions, windows):
    """ bCD-ring-integrity  - requires mask_list of 4-dim """
    r = restraints.DAT_restraint()
    r.topology = structure
    r.auto_apr = True
    r.continuous_apr = True
    r.amber_index = False
        
    r.mask1 = mask_list[0]
    r.mask2 = mask_list[1]
    r.mask3 = mask_list[2]
    r.mask4 = mask_list[3]
    f = lambda x: 108 if 'O5' in x else -114
    r.attach["target"] = f(mask_list[0])            # Degrees
    r.attach["fraction_list"] = attach_fractions
    r.attach["fc_final"] = 100.0                    # kcal/mol/radian**2
    r.pull["target_final"] = f(mask_list[0])        # Degrees
    r.pull["num_windows"] = windows[1]
    
    r.release['target'] = f(mask_list[0])   
    r.release['fraction_list'] = release_fractions
    r.release['fc_final'] = 0
    
    r.initialize()
    return r
  
  def build_host_conf_mask():
    conf_Mask = list()
    for i in range(1, 8):
        f = lambda x: 1 if x==8 else x
        conf_Mask.append([f":{f(i+1)},@O5", f":{f(i+1)},@C1", f":{i},@O4", f":{i}@C4"])
        conf_Mask.append([f":{f(i+1)},@C1", f":{i},@O4", f":{i},@C4", f":{i}@C5"])
    return conf_Mask
  
