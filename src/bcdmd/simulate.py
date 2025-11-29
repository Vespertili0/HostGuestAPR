import os, datetime, time, logging
import xml.etree.ElementTree as ET
from openmm.app import *
from openmm import *
from simtk.unit import *

class SimulationMD:
    """ Class to handle MD simulations for APR windows """
    def __init__(self, work_dir, window, step_size=2):
        self.WD = work_dir
        self.window = window
        self.folder = os.path.join(work_dir, window)
        self.step_size = step_size
        
    def _initialize_system(self, xml_file='system.xml', pdb_file='system.pdb'):
        self.FILES = os.listdir(self.folder)        
        with open(os.path.join(self.folder, xml_file), 'r') as file:
            self.system = XmlSerializer.deserialize(file.read())
        self.coords = app.PDBFile(os.path.join(self.folder, pdb_file))
        self.integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, self.step_size*femtoseconds)
        
    def _build_simulation(self, state):
        self.simulation = app.Simulation(self.coords.topology, self.system, self.integrator, state=state, 
                                         platform=Platform.getPlatformByName('CUDA'),
                                         platformProperties={'Precision': 'mixed'})
        if state is None:
            self.simulation.context.setPositions(self.coords.positions)
        print('simulation built')
            
    def _write_state(self, name):
        self.simulation.saveState(os.path.join(self.folder, f'{name}_state.txt'))
        
    def _prepare_MD(self, state, reportInterval, name, write_traj=True, restart=False):
        self._build_simulation(state=state)
        self.simulation.reporters.append(CheckpointReporter(os.path.join(self.folder, 'checkpnt.chk'), 500))
        if restart is True:
            DataFile = open(os.path.join(self.folder, f'{name}.txt'), 'a')
        if restart is False:
            DataFile = os.path.join(self.folder, f'{name}.txt')
        self.simulation.reporters.append(StateDataReporter(file=DataFile, reportInterval=reportInterval, 
                                                           step=True, time=True, potentialEnergy=True,
                                                           kineticEnergy=True, totalEnergy=True,
                                                           temperature=True, density=True, volume=True))
        if write_traj is True:
            self.simulation.reporters.append(DCDReporter(file=os.path.join(self.folder, f'{name}.dcd'), append=restart, 
                                                         reportInterval=reportInterval, enforcePeriodicBox=True))
        print('MD prepared')
    
    def minimize_anneal(self, target_T):
        """ Minimize and anneal the system to target temperature """
        self._initialize_system()
        self._build_simulation(state=None)
        self.simulation.minimizeEnergy(tolerance=5.0*kilojoules_per_mole, maxIterations=20000)
        annealing_increments = int(0.5 / self.step_size * 1e6 / 100) # 0.1 ns of 100 equal increments
        self.system.addForce(MonteCarloBarostat(1*bar, target_T*kelvin))
        T_factor = target_T / 100
        for i in range(100):
            self.integrator.setTemperature(T_factor * (100-i)*kelvin)
            self.simulation.step(annealing_increments)
        self._write_state(name='MIN')
        
    def runMD(self, eq_prod, ensemble, run_ns=1.5, write_traj=True, restart=False):
        """ Run MD simulation for equilibration or production """
        if self.step_size == 4:  ### write frames every 2ps ###
            snapshot = 500
        if self.step_size == 2:
            snapshot = 1000
        run_steps = int(run_ns / self.step_size * 1e6)
        
        self._initialize_system()
        if ensemble == 'NPT':
            self.system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
        if eq_prod == 'eq':
            state = os.path.join(self.folder, f'MIN_state.txt')
        if eq_prod == 'prod':
            state = os.path.join(self.folder, f'{ensemble}_eq_state.txt')
        file_name = f'{ensemble}_{eq_prod}'
        no_attempt = 0
        if restart is True:
            try:
                self._prepare_MD(state=state, reportInterval=snapshot, name=file_name,
                                 restart=True, write_traj=write_traj)
                self.simulation.loadCheckpoint(os.path.join(self.folder, 'checkpnt.chk'))
                print(f'chk loaded ... continue')
                self.simulation.step(run_steps)
            except OpenMMException:
                print(f'...{self.window} failed restart')
        if restart is False:
            self._prepare_MD(state=state, reportInterval=snapshot, name=file_name,
                 write_traj=write_traj)
            self.simulation.step(run_steps)
        self._write_state(name=file_name)


def MDrunAPR(SIM_LIST, work_dir):
    """ Run MD simulations for all APR windows """
    FAILED = list()
    for i, window in enumerate(SIM_LIST):
        try:
            if 'p' in window:
                prod_ns = 3   # normal is 10ns
                if i > 0:
                    print('cooldown')
                    time.sleep(60)
            else:
                prod_ns = 2
                if i > 0 and i % 3 == 0:
                    print('cooldown')
                    time.sleep(30)
            print(f'...{window}...')      
            apr = SimulationMD(work_dir=work_dir, window=window, step_size=4)
            apr.minimize_anneal(300)
            print(datetime.datetime.now(), '...minimize & annealing done')
            apr.runMD(eq_prod='eq', ensemble='NPT', write_traj=False)
            apr.runMD(eq_prod='prod', ensemble='NPT', run_ns=prod_ns, write_traj=True, restart=False)
            print(datetime.datetime.now(), '...MD simulation done')
        except:
            FAILED.append(window)
            print(datetime.datetime.now(), f'\n...{window} failed...\n')
    return FAILED
