import os
import numpy as np
import pandas as pd
import parmed as pmd
import paprika.analysis as analysis
from paprika.io import load_restraints
from paprika.restraints.utils import extract_guest_restraints

class APRanalysis:
  def __init__(self, work_dir, guest_resname):
    system_restraints = load_restraints(f'{work_dir}/restraints.json')
    self.free_energy = analysis.fe_calc()
    self.free_energy.prmtop = "system.pdb"
    self.free_energy.trajectory = "NPT_prod.dcd"
    self.free_energy.path = f'{work_dir}/apr_windows/'
    self.free_energy.restraint_list = self.system_restraints
    self.free_energy.collect_data()
    
    structure = pmd.load_file(f'{work_dir}/complex/complex.prmtop', f'{work_dir}/complex/complex.rst7', structure=True)
    self.guest_restraints = extract_guest_restraints(structure, system_restraints, guest_resname=guest_resname)
    self.free_energy.compute_ref_state_work(self.guest_restraints)
    
  def analyze(self, temperature=300, ti_bootcycles=100000):
    self.free_energy.methods = ['mbar-autoc', 'mbar-block', 'ti-block']
    self.free_energy.temperature = temperature
    self.free_energy.ti_matrix = "diagonal"
    self.free_energy.bootcycles = ti_bootcycles
    self.free_energy.compute_free_energy()
  
  def summarize(self):
    DATA = pd.DataFrame(columns=['method', 'dG_total', 'SEM_total', 'attach_G', 'attach_SEM',
                                 'pull_G', 'pull_SEM ', 'release_G', 'release_SEM'])
    for method in self.free_energy.methods:
        a_G = self.free_energy.results["attach"][method]["fe"]
        a_sem = self.free_energy.results["attach"][method]["sem"]
        p_G = self.free_energy.results["pull"][method]["fe"]
        p_sem = self.free_energy.results["pull"][method]["sem"]
        try:
          r_G = self.free_energy.results["release"][method]["fe"]
          r_sem = self.free_energy.results['release'][method]['sem']
        except:
          r_G, r_sem = 0, 0
        ref_G = self.free_energy.results["ref_state_work"]
        
        binding_affinity = -1 * ( a_G + p_G - r_G + ref_G)
        sem = np.sqrt( a_sem**2 + p_sem**2 + r_sem**2)
        new = pd.DataFrame([[method, binding_affinity, sem, a_G, a_sem,
                            p_G, p_sem, r_G, r_sem]], columns=DATA.columns)
        DATA = pd.concat([DATA, new], ignore_index=True)
    return DATA
