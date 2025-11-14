# Figure - Simulation of CHO fedbatch reactor with well designed feed profile
#          with functions added to facilitate explorative simulation work 
#
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2020-02-24 - Python3 script for Windows and Linux
#            - Added system_info() that prints system information
#            - Change newplot and simu using objectoriented diagrams 
#            - Simplified handling of simulation results
#            - Tested with JModelica 2.14 and seems ok
#            - Introduced locale and setting of it - important for OpenModelica
#            - Introduced check of platform to adapt code for Windows/Linux
#            - Correct print() and np.nan
#            - Adpated for Jupyter
# 2020-02-28 - Corrected parSet() VX_0 should be VXv_0
#------------------------------------------------------------------------------------------------------------------
# 2020-06-18 - Adapted for BR5m
#------------------------------------------------------------------------------------------------------------------
# 2020-07-15 - Adapted for BP6a
# 2020-07-18 - Eliminated use of get() before simulation to comply with FMU standard
# 2020-07-27 - Introduce choice of Linux FMU - JModelica or OpenModelica
# 2020-07-28 - Change of simu('cont') and handling of stateDict and model.get..
# 2020-11-21 - Adapted to ReactorType with n_inlets, n_outlets and n_ports
# 2021-02-05 - Adjust describe() for change to liquidphase
#------------------------------------------------------------------------------------------------------------------
# 2021-02-11 - Adapt for coming BPL_v2
# 2021-02-16 - Adapt for further restructing in packages and later divide into files
# 2021-03-20 - Adapt for BPL ver 2.0.3
# 2021-04-15 - Adapt for Fedbatch2 with use of MSL CombiTimaTable
# 2021-04-23 - Adapt for BPL ver 2.0.4
# 2021-07-02 - Modify interaction to the current state - now application part small and general functions ok
# 2021-08-05 - Introduced describe_parts() and corrected disp() to handle number of displayed decimals 
# 2021-09-13 - Tested with BPL ver 2.0.7
# 2021-10-02 - Updated system_info() with FMU-explore version
# 2021-11-21 - Included bioreactor.broth_decay also after change of name
# 2022-01-26 - Updated to FMU-explore 0.8.8
# 2022-02-01 - Updated to FMU-explore 0.8.9
# 2022-03-28 - Updated to FMU-explore 0.9.0 - model.reset(), and par(), init()
# 2022-08-19 - Updated for BPL ver 2.1.0 beta and FMU-exolre 0.9.2
# 2022-10-06 - Updated for FMU-explore 0.9.5 with disp() that do not include extra parameters with parLocation
# 2023-02-15 - Consolidate FMU-explore to 0.9.6 and means parCheck and par() udpate and simu() with opts as arg
# 2023-03-29 - Update FMU-explore 0.9.7
# 2024-03-01 - Update FMU-explore 0.9.8 and 0.9.9 - now with _0 replaced with _start everywhere
# 2024-10-16 - Update FMU-explore 1.0.0
# 2024-11-04 - Removed GUI from process diagram
# 2024-11-07 - Update BPL 2.3.0
# 2025-02-20 - Try with CHO - extended with Xl i.e. lysed cells that bring toxicity
# 2025-04-15 - Added describe('process') using model.get_description() and the show text string for the main model
# 2025-06-12 - Test MSL 4.1.0 with OpenModelica genreated FMU
# 2025-07-23 - Update BPL 2.3.1
# 2025-11-08 - FMU-explore 1.0.2
# 2025-11-13 - Test FMU-explore 1.0.1h and global declaration removed outside functions
# 2025-11-14 - FMU-explore 1.0.2 corrected
#-------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.image as img
import zipfile 
 
from pyfmi import load_fmu
from pyfmi.fmi import FMUException

from itertools import cycle
from importlib_metadata import version   # included in future Python 3.8

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   flag_vendor = 'JM'
   flag_type = 'CS'
   fmu_model ='BPL_CHO_Fedbatch_windows_jm_cs.fmu'        
   model = load_fmu(fmu_model, log_level=0)  
elif platform.system() == 'Linux':   
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-compiled OpenModelica') 
      if flag_type in ['CS','cs']:         
         fmu_model ='BPL_CHO_Fedbatch_linux_om_cs.fmu'    
         model = load_fmu(fmu_model, log_level=0) 
      if flag_type in ['ME','me']:         
         fmu_model ='BPL_CHO_Fedbatch_linux_om_me.fmu'    
         model = load_fmu(fmu_model, log_level=0)
   else:    
      print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = model.simulate_options()
   opts_std['silent_mode'] = True
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary'     
elif flag_type in ['ME', 'me']:
   opts_std = model.simulate_options()
   opts_std["CVode_options"]["verbosity"] = 50 
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary'  
else:    
   print('There is no FMU for this platform')
  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
   MSL_usage = model.get('MSL.usage')[0]
   MSL_version = model.get('MSL.version')[0]
   BPL_version = model.get('BPL.version')[0]
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '4.1.0 - used components: RealInput, RealOutput. CombiTimeTable, Types' 
   MSL_version = '4.1.0'
   BPL_version = 'Bioprocess Library version 2.3.1' 
else:    
   print('There is no FMU for this platform')

# Simulation time
simulationTime = 120.0
prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture', 'bioreactor.broth_decay']

# Process diagram
fmu_process_diagram ='BPL_CHO_Fedbatch_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, newplot(), describe()
#------------------------------------------------------------------------------------------------------------------
   
# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = model.get_states_list()
stateValue.update(timeDiscreteStates)

# Create parValue
parValue = {}
parValue['V_start']    = 0.35          # L
parValue['VXv_start'] = 0.35*0.2       
parValue['VXd_start'] = 0.0  
parValue['VXl_start'] = 0.0              
parValue['VG_start'] = 0.35*18.0       
parValue['VGn_start'] = 0.35*2.4       
parValue['VL_start'] = 0.0             
parValue['VN_start'] = 0.0             

parValue['qG_max1'] = 0.2971
parValue['qG_max2'] = 0.0384
parValue['qGn_max1'] = 0.1238
parValue['qGn_max2'] = 0.0218
parValue['mu_d_max'] = 0.1302
parValue['k_toxic'] = 0.0
parValue['alpha'] = 0
parValue['beta'] = 10.0/24

parValue['k_lysis_v'] = 0.0
parValue['k_lysis_d'] = 0.0

parValue['G_in']  =  15.0          # mM
parValue['Gn_in']  =  4.0          # mM
parValue['t0'] =   0.0             # h
parValue['F0'] =   0.0             # L/h
parValue['t1'] =  50.0             # h
parValue['F1'] =   0.0012          # L/h
parValue['t2'] =  64.0             # h
parValue['F2'] =   0.0020          # L/h
parValue['t3'] =  78.0             # h
parValue['F3'] =   0.0040          # L/h
parValue['t4'] =  92.0             # h
parValue['F4'] =   0.0080          # L/h
parValue['t5'] = 106.0             # h
parValue['F5'] =   0.012           # L/h
parValue['t6'] = 150.0             # h
parValue['F6'] =   0.012           # L/h

parLocation = {}
parLocation['V_start'] = 'bioreactor.V_start'
parLocation['VXv_start'] = 'bioreactor.m_start[1]'
parLocation['VXd_start'] = 'bioreactor.m_start[2]'
parLocation['VXl_start'] = 'bioreactor.m_start[3]'
parLocation['VG_start'] = 'bioreactor.m_start[4]'
parLocation['VGn_start'] = 'bioreactor.m_start[5]'
parLocation['VL_start'] = 'bioreactor.m_start[6]'
parLocation['VN_start'] = 'bioreactor.m_start[7]'

parLocation['qG_max1'] = 'bioreactor.culture.qG_max1'
parLocation['qG_max2'] = 'bioreactor.culture.qG_max2'
parLocation['qGn_max1'] = 'bioreactor.culture.qGn_max1'
parLocation['qGn_max2'] = 'bioreactor.culture.qGn_max2'
parLocation['mu_d_max'] = 'bioreactor.culture.mu_d_max'
parLocation['k_toxic'] = 'bioreactor.culture.k_toxic'
parLocation['alpha'] = 'bioreactor.culture.alpha'
parLocation['beta'] = 'bioreactor.culture.beta'

parLocation['k_lysis_v'] = 'bioreactor.broth_decay.k_lysis_v'
parLocation['k_lysis_d'] = 'bioreactor.broth_decay.k_lysis_d'

parLocation['G_in'] = 'feedtank.c_in[4]'
parLocation['Gn_in'] = 'feedtank.c_in[5]'

parLocation['t0'] = 'dosagescheme.table[1,1]'
parLocation['F0'] = 'dosagescheme.table[1,2]'
parLocation['t1'] = 'dosagescheme.table[2,1]'
parLocation['F1'] = 'dosagescheme.table[2,2]'
parLocation['t2'] = 'dosagescheme.table[3,1]'
parLocation['F2'] = 'dosagescheme.table[3,2]'
parLocation['t3'] = 'dosagescheme.table[4,1]'
parLocation['F3'] = 'dosagescheme.table[4,2]'
parLocation['t4'] = 'dosagescheme.table[5,1]'
parLocation['F4'] = 'dosagescheme.table[5,2]'
parLocation['t5'] = 'dosagescheme.table[6,1]'
parLocation['F5'] = 'dosagescheme.table[6,2]'
parLocation['t6'] = 'dosagescheme.table[7,1]'
parLocation['F6'] = 'dosagescheme.table[7,2]'

# Extra only for describe()
parLocation['mu'] = 'bioreactor.culture.mu'
parLocation['mu_d'] = 'bioreactor.culture.mu_d'

# Parameter value check - especially for hysteresis to avoid runtime error
parCheck = []
parCheck.append("parValue['V_start'] > 0")
parCheck.append("parValue['VXv_start'] >= 0")
parCheck.append("parValue['VG_start'] >= 0")
parCheck.append("parValue['VGn_start'] >= 0")
parCheck.append("parValue['VL_start'] >= 0")
parCheck.append("parValue['VN_start'] >= 0")

# Create list of diagrams to be plotted by simu()
diagrams = []

def newplot(title='Fedbatch cultivation',  plotType='TimeSeries'):
   """ Standard plot window,
        title = '' """
   
   # Reset pens
   setLines()
   
   # Globals
   global ax11, ax12, ax13, ax21, ax22, ax23, ax31, ax32, ax33, ax41, ax42, ax43, ax51    
    
   # Plot diagram 
   if plotType == 'TimeSeries':

      plt.figure()
      ax11 = plt.subplot(4,2,1); ax12 = plt.subplot(4,2,2)
      ax21 = plt.subplot(4,2,3); ax22 = plt.subplot(4,2,4)
      ax31 = plt.subplot(4,2,5); ax32 = plt.subplot(4,2,6)
      ax41 = plt.subplot(4,2,7); ax42 = plt.subplot(4,2,8)    

      ax11.set_title(title)
      ax11.grid()
      ax11.set_ylabel('Glucose conc [mM]')

      ax12.grid()
      ax12.set_ylabel('Lactate conc [mM]')

      ax21.grid()
      ax21.set_ylabel('Glutamine conc [mM]')

      ax22.grid()
      ax22.set_ylabel('Ammonia conc [mM]')

      ax31.grid()
      ax31.set_ylabel('Viable cells conc [1E6/mL]')

      ax32.grid()
      ax32.set_ylabel('Dead cells conc [1E6/mL]')

      ax41.grid()
      ax41.set_ylabel('Feed rate [L/h]')
      ax41.set_xlabel('Time [h]')

      ax42.grid()
      ax42.set_ylabel('Volume [L]')
      ax42.set_xlabel('Time [h]')

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(t,sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")        
      diagrams.append("ax41.plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax42.plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 
   
   if plotType == 'TimeSeries1':

      plt.figure()
      ax11 = plt.subplot(4,2,1); ax12 = plt.subplot(4,2,2)
      ax21 = plt.subplot(4,2,3); ax22 = plt.subplot(4,2,4)
      ax31 = plt.subplot(4,2,5); ax32 = plt.subplot(4,2,6)
      ax41 = plt.subplot(4,2,7); ax42 = plt.subplot(4,2,8)    

      ax11.set_title(title)
      ax11.grid()
      ax11.set_ylabel('Glucose conc [mM]')

      ax12.grid()
      ax12.set_ylabel('Lactate conc [mM]')

      ax21.grid()
      ax21.set_ylabel('Glutamine conc [mM]')

      ax22.grid()
      ax22.set_ylabel('Ammonia conc [mM]')

      ax31.grid()
      ax31.set_ylabel('Viable cells conc [1E6/mL]')

      ax32.grid()
      ax32.set_ylabel('Dead and lysed cells [1E6/mL]')

      ax41.grid()
      ax41.set_ylabel('Feed rate [L/h]')
      ax41.set_xlabel('Time [h]')

      ax42.grid()
      ax42.set_ylabel('Volume [L]')
      ax42.set_xlabel('Time [h]')

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(t,sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")  
      diagrams.append("ax32.plot(t,sim_res['bioreactor.c[3]'], color='k', linestyle=linetype)")       
      diagrams.append("ax41.plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax42.plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

    
   # Plot diagram 
   if plotType == 'TimeSeries2':
   
      plt.figure()
      ax11 = plt.subplot(4,2,1); ax12 = plt.subplot(4,2,2)
      ax21 = plt.subplot(4,2,3); ax22 = plt.subplot(4,2,4)
      ax31 = plt.subplot(4,2,5); ax32 = plt.subplot(4,2,6)
      ax41 = plt.subplot(4,2,7); ax42 = plt.subplot(4,2,8)    

      ax11.set_title(title)
      ax11.grid()
      ax11.set_ylabel('Glucose conc [mM]')

      ax12.grid()
      ax12.set_ylabel('Lactate conc [mM]')

      ax21.grid()
      ax21.set_ylabel('Glutamine conc [mM]')

      ax22.grid()
      ax22.set_ylabel('Ammonia conc [mM]')

      ax31.grid()
      ax31.set_ylabel('Viable cells [1E6]')

      ax32.grid()
      ax32.set_ylabel('Dead and lysed cells [1E6/mL]')

      ax41.grid()
      ax41.set_ylabel('Feed rate [L/h]')
      ax41.set_xlabel('Time [h]')

      ax42.grid()
      ax42.set_ylabel('Volume [L]')
      ax42.set_xlabel('Time [h]')

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(t,sim_res['bioreactor.m[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")  
      diagrams.append("ax32.plot(t,sim_res['bioreactor.c[3]'], color='k', linestyle=linetype)")       
      diagrams.append("ax41.plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax42.plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 
      
   if plotType == 'Textbook_3':
 
      plt.figure()
      ax11 = plt.subplot(5,3,1); ax12 = plt.subplot(5,3,2); ax13 = plt.subplot(5,3,3)
      ax21 = plt.subplot(5,3,4); ax22 = plt.subplot(5,3,5); ax23 = plt.subplot(5,3,6)
      ax31 = plt.subplot(5,3,7); ax32 = plt.subplot(5,3,8); ax33 = plt.subplot(5,3,9)
      ax41 = plt.subplot(5,3,10); ax42 = plt.subplot(5,3,11); ax43 = plt.subplot(5,3,12)
      ax51 = plt.subplot(5,3,13) 

      ax11.set_title(title)
      ax11.grid()
      ax11.set_ylabel('Glucose [mM]')

      ax12.grid()
      ax12.set_ylabel('Lactate [mM]')
      
      ax13.grid()
      ax13.set_ylabel('qG []')

      ax21.grid()
      ax21.set_ylabel('Glutamine [mM]')

      ax22.grid()
      ax22.set_ylabel('Ammonia [mM]')
      
      ax23.grid()
      ax23.set_ylabel('qGn []')

      ax31.grid()
      ax31.set_ylabel('Viable cell [1E6/mL]')

      ax32.grid()
      ax32.set_ylabel('Dead cell conc [1E6/mL]')
      ax32.set_xlabel('Time [h]')
      
      ax33.grid()
      ax33.set_ylabel('mu [1/h]')
      ax33.set_xlabel('Time [h]')

      ax41.grid()
      ax41.set_ylabel('Feed rate [L/h]')

      ax42.grid()
      ax42.set_ylabel('mAb []')

      ax43.grid()
      ax43.set_ylabel('qP []')

      ax51.grid()
      ax51.set_ylabel('Volume [L]')      
      ax51.set_xlabel('Time [h]')
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(t,sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")       
      diagrams.append("ax41.plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)") 
      diagrams.append("ax42.plot(t,sim_res['bioreactor.c[8]'], color='g', linestyle=linetype)")       
      diagrams.append("ax51.plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

      diagrams.append("ax13.set_title('- cell specific rates')") 
      diagrams.append("ax13.plot(t,-(sim_res['bioreactor.culture.q[4]']+sim_res['bioreactor.culture.qG_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax13.plot(t,-sim_res['bioreactor.culture.q[4]'], color='b', linestyle=linetype)") 
      diagrams.append("ax23.plot(t,-(sim_res['bioreactor.culture.q[5]']+sim_res['bioreactor.culture.qGn_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax23.plot(t,-sim_res['bioreactor.culture.q[5]'], color='b', linestyle=linetype)") 
      diagrams.append("ax33.plot(t,sim_res['bioreactor.culture.q[1]'], color='b', linestyle=linetype)") 
      diagrams.append("ax43.plot(t,sim_res['bioreactor.culture.q[8]'], color='g', linestyle=linetype)") 
      
      diagrams.append("ax11.set_ylim(0)")
      diagrams.append("ax13.set_ylim(0)")
      diagrams.append("ax32.set_ylim(ax31.get_ylim())")      

      
def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""

   if name == 'process':
      model.get_description()

   elif name == 'culture':
      print('Reactor culture CHO-MAb - cell line HB-58 American Culture Collection ATCC') 

   elif name in ['broth', 'liquidphase', 'liquid-phase']:

      Xv  = model.get('liquidphase.Xv')[0]; 
      Xv_description = model.get_variable_description('liquidphase.Xv'); 
      Xv_mw = model.get('liquidphase.mw[1]')[0]
      
      Xd = model.get('liquidphase.Xd')[0]; 
      Xd_description = model.get_variable_description('liquidphase.Xd'); 
      Xd_mw = model.get('liquidphase.mw[2]')[0]

      Xl = model.get('liquidphase.Xl')[0]; 
      Xl_description = model.get_variable_description('liquidphase.Xl'); 
      Xl_mw = model.get('liquidphase.mw[3]')[0]
      
      G = model.get('liquidphase.G')[0]; 
      G_description = model.get_variable_description('liquidphase.G'); 
      G_mw = model.get('liquidphase.mw[4]')[0]
      
      Gn = model.get('liquidphase.Gn')[0]; 
      Gn_description = model.get_variable_description('liquidphase.Gn'); 
      Gn_mw = model.get('liquidphase.mw[5]')[0]
      
      L = model.get('liquidphase.L')[0]; 
      L_description = model.get_variable_description('liquidphase.L'); 
      L_mw = model.get('liquidphase.mw[6]')[0]
      
      N = model.get('liquidphase.N')[0]; 
      N_description = model.get_variable_description('liquidphase.N'); 
      N_mw = model.get('liquidphase.mw[7]')[0]
      
      Pr = model.get('liquidphase.Pr')[0]; 
      Pr_description = model.get_variable_description('liquidphase.Pr'); 
      Pr_mw = model.get('liquidphase.mw[8]')[0]

      print('Reactor broth substances included in the model')
      print()
      print(Xv_description, 'index = ', Xv, 'molecular weight = ', Xv_mw, 'Da')
      print(Xd_description, '  index = ', Xd, 'molecular weight = ', Xd_mw, 'Da')
      print(Xl_description, ' index = ', Xl, 'molecular weight = ', Xl_mw, 'Da')      
      print(G_description, '     index = ', G, 'molecular weight = ', G_mw, 'Da')
      print(Gn_description, '   index = ', Gn, 'molecular weight = ', Gn_mw, 'Da')
      print(L_description, '     index = ', L, 'molecular weight = ', L_mw, 'Da')
      print(N_description, '     index = ', N, 'molecular weight = ', N_mw, 'Da')
      print(Pr_description, '     index = ', Pr, 'molecular weight = ', Pr_mw, 'Da')

   elif name in ['parts']:
      describe_parts(component_list_minimum)

   elif name in ['MSL']:
      describe_MSL()

   else:
      describe_general(name, decimals)

#------------------------------------------------------------------------------------------------------------------
#  General code 
FMU_explore = 'FMU-explore version 1.0.2'
#------------------------------------------------------------------------------------------------------------------

# Define function par() for parameter update
def par(*x, parValue=parValue, **x_kwarg):
   """ Set parameter values if available in the predefined dictionaryt parValue. """
   x_kwarg.update(*x)
   x_temp = {}
   for key in x_kwarg.keys():
      if key in parValue.keys():
         x_temp.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an accessible parameter - check the spelling')
   parValue.update(x_temp)
   
   parErrors = [requirement for requirement in parCheck if not(eval(requirement))]
   if not parErrors == []:
      print('Error - the following requirements do not hold:')
      for index, item in enumerate(parErrors): print(item)

# Define function init() for initial values update
def init(*x, parValue=parValue, **x_kwarg):
   """ Set initial values and the name should contain string '_start' to be accepted.
       The function can handle general parameter string location names if entered as a dictionary. """
   x_kwarg.update(*x)
   x_init={}
   for key in x_kwarg.keys():
      if '_start' in key: 
         x_init.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an initial value, use par() instead - check the spelling')
   parValue.update(x_init)

# Define how to read dictionary for parameter values
def readParValue(file, sheet, parValue=parValue):
   """ Read parameter short names and values from an Excel-file from defined sheet. For use in the notebook!
       Return a dictionary."""
   parValue_local = {} 
   table = pd.ExcelFile(file).parse(sheet)
   for k in list(range(len(table))):
      parValue_local[table['Par'][k]] = table['Value'][k]
   parValue.update(parValue_local)

# Define how to read dictionary for parameter location
def readParLocation(file, parLocation=parLocation):
   """ Read parameter short and long names from an Excel-file sheet by sheet. For use in the notebook!
       Return a dictionary."""
   sheets = ['initial_values','feed_AB', 'feed_G', 'culture', 'broth_decay']
   parLocation_local = {}
   for sheet in sheets:
      table = pd.ExcelFile(file).parse(sheet)
      for k in list(range(len(table))):
         parLocation_local[table['Par'][k]] = table['Location'][k]
   parLocation.update(parLocation_local)
      
def disp(name='', decimals=3, mode='short', parValue=parValue, parLocation=parLocation):
   """ Display intial values and parameters in the model that include "name" and is in parLocation list.
       Note, it does not take the value from the dictionary par but from the model. """
   global model

   def dict_reverser(d):
      seen = set()
      return {v: k for k, v in d.items() if v not in seen or seen.add(v)}
   
   if mode in ['short']:
      k = 0
      for Location in [parLocation[k] for k in parValue.keys()]:
         if name in Location:
            if type(model.get(Location)[0]) != np.bool_:
               print(dict_reverser(parLocation)[Location] , ':', np.round(model.get(Location)[0],decimals))
            else:
               print(dict_reverser(parLocation)[Location] , ':', model.get(Location)[0])               
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parValue.keys():
            if name in parName:
               if type(model.get(Location)[0]) != np.bool_:
                  print(parName,':', np.round(model.get(parLocation[parName])[0],decimals))
               else: 
                  print(parName,':', model.get(parLocation[parName])[0])
   if mode in ['long','location']:
      k = 0
      for Location in [parLocation[k] for k in parValue.keys()]:
         if name in Location:
            if type(model.get(Location)[0]) != np.bool_:       
               print(Location,':', dict_reverser(parLocation)[Location] , ':', np.round(model.get(Location)[0],decimals))
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parValue.keys():
            if name in parName:
               if type(model.get(Location)[0]) != np.bool_:
                  print(parLocation[parName], ':', dict_reverser(parLocation)[Location], ':', parName,':', 
                     np.round(model.get(parLocation[parName])[0],decimals))

# Line types
def setLines(lines=['-','--',':','-.']):
   """Set list of linetypes used in plots"""
   global linecycler
   linecycler = cycle(lines)

# Show plots from sim_res, just that
def show(diagrams=diagrams):
   """Show diagrams chosen by newplot()"""
   # Plot pen
   linetype = next(linecycler)    
   # Plot diagrams 
   for command in diagrams: eval(command)

# Simulation
def simu(simulationTimeLocal=simulationTime, mode='Initial', options=opts_std, \
         diagrams=diagrams,timeDiscreteStates=timeDiscreteStates, stateValue=stateValue, \
         parValue=parValue, parLocation=parLocation, fmu_model=fmu_model):         
   """Model loaded and given intial values and parameter before,
      and plot window also setup before."""
    
   # Global variables
   global model, prevFinalTime, sim_res, t
   
   # Simulation flag
   simulationDone = False
   
   # Transfer of argument to global variable
   simulationTime = simulationTimeLocal 
      
   # Check parValue
   value_missing = 0
   for key in parValue.keys():
      if parValue[key] in [np.nan, None, '']:
         print('Value missing:', key)
         value_missing =+1
   if value_missing>0: return
         
   # Load model
   if model is None:
      model = load_fmu(fmu_model) 
   model.reset()
      
   # Run simulation
   if mode in ['Initial', 'initial', 'init']:
      # Set parameters and intial state values:
      for key in parValue.keys():
         model.set(parLocation[key],parValue[key])   
      # Simulate
      sim_res = model.simulate(final_time=simulationTime, options=options)  
      simulationDone = True
   elif mode in ['Continued', 'continued', 'cont']:

      if prevFinalTime == 0: 
         print("Error: Simulation is first done with default mode = init'")      
      else:
         
         # Set parameters and intial state values:
         for key in parValue.keys():
            model.set(parLocation[key],parValue[key])                

         for key in stateValue.keys():
            if not key[-1] == ']':
               if key[-3:] == 'I.y': 
                  model.set(key[:-10]+'I_start', stateValue[key]) 
               elif key[-3:] == 'D.x': 
                  model.set(key[:-10]+'D_start', stateValue[key]) 
               else:
                  model.set(key+'_start', stateValue[key])
            elif key[-3] == '[':
               model.set(key[:-3]+'_start'+key[-3:], stateValue[key]) 
            elif key[-4] == '[':
               model.set(key[:-4]+'_start'+key[-4:], stateValue[key]) 
            elif key[-5] == '[':
               model.set(key[:-5]+'_start'+key[-5:], stateValue[key]) 
            else:
               print('The state vecotr has more than 1000 states')
               break

         # Simulate
         sim_res = model.simulate(start_time=prevFinalTime,
                                 final_time=prevFinalTime + simulationTime,
                                 options=options) 
         simulationDone = True             
   else:
      print("Simulation mode not correct")

   if simulationDone:
    
      # Extract data
      t = sim_res['time']
 
      # Plot diagrams
      linetype = next(linecycler)    
      for command in diagrams: eval(command)
            
      # Store final state values stateValue:
      for key in list(stateValue.keys()): stateValue[key] = model.get(key)[0]        

      # Store time from where simulation will start next time
      prevFinalTime = model.time
   
   else:
      print('Error: No simulation done')
      
# Describe model parts of the combined system
def describe_parts(component_list=[]):
   """List all parts of the model""" 
       
   def model_component(variable_name):
      i = 0
      name = ''
      finished = False
      if not variable_name[0] == '_':
         while not finished:
            name = name + variable_name[i]
            if i == len(variable_name)-1:
                finished = True 
            elif variable_name[i+1] in ['.', '(']: 
                finished = True
            else: 
                i=i+1
      if name in ['der', 'temp_1', 'temp_2', 'temp_3', 'temp_4', 'temp_5', 'temp_6', 'temp_7']: name = ''
      return name
    
   variables = list(model.get_model_variables().keys())
        
   for i in range(len(variables)):
      component = model_component(variables[i])
      if (component not in component_list) \
      & (component not in ['','BPL', 'Customer', 'today[1]', 'today[2]', 'today[3]', 'temp_2', 'temp_3']):
         component_list.append(component)
      
   print(sorted(component_list, key=str.casefold))
   
def describe_MSL(flag_vendor=flag_vendor):
   """List MSL version and components used"""
   print('MSL:', MSL_usage)
 
# Describe parameters and variables in the Modelica code
def describe_general(name, decimals, parLocation=parLocation):
  
   if name == 'time':
      description = 'Time'
      unit = 'h'
      print(description,'[',unit,']')
      
   elif name in parLocation.keys():
      description = model.get_variable_description(parLocation[name])
      value = model.get(parLocation[name])[0]
      try:
         unit = model.get_variable_unit(parLocation[name])
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)            
      else:
        print(description, ':', np.round(value, decimals), '[',unit,']')
                  
   else:
      description = model.get_variable_description(name)
      value = model.get(name)[0]
      try:
         unit = model.get_variable_unit(name)
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)     
      else:
         print(description, ':', np.round(value, decimals), '[',unit,']')
         
# Plot process diagram
def process_diagram(fmu_model=fmu_model, fmu_process_diagram=fmu_process_diagram):   
   try:
       process_diagram = zipfile.ZipFile(fmu_model, 'r').open('documentation/processDiagram.png')
   except KeyError:
       print('No processDiagram.png file in the FMU, but try the file on disk.')
       process_diagram = fmu_process_diagram
   try:
       plt.imshow(img.imread(process_diagram))
       plt.axis('off')
       plt.show()
   except FileNotFoundError:
       print('And no such file on disk either')
         
# Describe framework
def BPL_info():
   print()
   print('Model for the process has been setup. Key commands:')
   print(' - par()       - change of parameters and initial values')
   print(' - init()      - change initial values only')
   print(' - simu()      - simulate and plot')
   print(' - newplot()   - make a new plot')
   print(' - show()      - show plot from previous simulation')
   print(' - disp()      - display parameters and initial values from the last simulation')
   print(' - describe()  - describe culture, broth, parameters, variables with values/units')
   print()
   print('Note that both disp() and describe() takes values from the last simulation')
   print('and the command process_diagram() brings up the main configuration')
   print()
   print('Brief information about a command by help(), eg help(simu)') 
   print('Key system information is listed with the command system_info()')

def system_info():
   """Print system information"""
   FMU_type = model.__class__.__name__
   print()
   print('System information')
   print(' -OS:', platform.system())
   print(' -Python:', platform.python_version())
   try:
       scipy_ver = scipy.__version__
       print(' -Scipy:',scipy_ver)
   except NameError:
       print(' -Scipy: not installed in the notebook')
   print(' -PyFMI:', version('pyfmi'))
   print(' -FMU by:', model.get_generation_tool())
   print(' -FMI:', model.get_version())
   print(' -Type:', FMU_type)
   print(' -Name:', model.get_name())
   print(' -Generated:', model.get_generation_date_and_time())
   print(' -MSL:', MSL_version)    
   print(' -Description:', BPL_version)   
   print(' -Interaction:', FMU_explore)
   
#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

BPL_info()