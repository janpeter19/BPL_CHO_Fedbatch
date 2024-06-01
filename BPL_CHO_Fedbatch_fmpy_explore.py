# Figure - Simulation of CHO fedbatch reactor with well designed feed profile
#          with functions added to facilitate explorative simulation work 
#
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2022-10-05 - Updated for FMU-explore 0.9.5 with disp() that do not include extra parameters with parLocation
# 2023-02-08 - Updated to FMU-explore 0.9.6e
# 2023-02-13 - Consolidate FMU-explore to 0.9.6 and means parCheck and par() udpate and simu() with opts as arg
# 2023-02-28 - Update FMU-explore for FMPy 0.9.6 in one leap and added list key_variables for logging
# 2023-03-22 - Correcting the script by including logging of states in a pedestrian way
# 2023-03-23 - Update FMU-explore to 0.9.7c
# 2023-03-28 - Update FMU-explore 0.9.7
# 2023-04-21 - Compiled for Ubuntu 20.04 and changed BPL_version
# 2023-05-31 - Adjusted to from importlib.meetadata import version
# 2023-09-12 - Updated to FMU-explore 0.9.8 and introduced process diagram
# 2024-03-07 - Update FMU-explore 0.9.9 - now with _0 replaced with _start everywhere
# 2024-05-20 - Updated the OpenModelica version to 1.23.0-dev
# 2024-06-01 - Corrected model_get() to handle string values as well - improvement very small and keep ver 1.0.0
#-------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.image as img
import zipfile 

from fmpy import simulate_fmu
from fmpy import read_model_description
import fmpy as fmpy

from itertools import cycle
from importlib.metadata import version  

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
global fmu_model, model
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   fmu_model ='BPL_CHO_Fedbatch_windows_jm_cs.fmu'        
   model_description = read_model_description(fmu_model)  
   flag_vendor = 'JM'
   flag_type = 'CS' 
elif platform.system() == 'Linux':  
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-comiled OpenModelica 1.23.0-dev') 
      if flag_type in ['CS','cs']:         
         fmu_model ='BPL_CHO_Fedbatch_linux_om_cs.fmu'    
         model_description = read_model_description(fmu_model) 
      if flag_type in ['ME','me']:         
         fmu_model ='BPL_CHO_Fedbatch_linux_om_me.fmu'    
         model_description = read_model_description(fmu_model) 
   else:    
      print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = {'ncp': 500}
elif flag_type in ['ME', 'me']:
   opts_std = {'ncp': 500}
else:    
   print('There is no FMU for this platform')

  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
   constants = [v for v in model_description.modelVariables if v.causality == 'local'] 
   MSL_usage = [x[1] for x in [(constants[k].name, constants[k].start) \
                     for k in range(len(constants))] if 'MSL.usage' in x[0]][0]   
   MSL_version = [x[1] for x in [(constants[k].name, constants[k].start) \
                       for k in range(len(constants))] if 'MSL.version' in x[0]][0]
   BPL_version = [x[1] for x in [(constants[k].name, constants[k].start) \
                       for k in range(len(constants))] if 'BPL.version' in x[0]][0] 
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '3.2.3 - used components: RealInput, RealOutput, CombiTimeTable, Types' 
   MSL_version = '3.2.3'
   BPL_version = 'Bioprocess Library version 2.2.0' 
else:    
   print('There is no FMU for this platform')
   
# Simulation time
global simulationTime; simulationTime = 120.0
global prevFinalTime; prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture', 'bioreactor.broth_decay']

# Provide process diagram on disk
fmu_process_diagram ='BPL_GUI_CHO_Fedbatch_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateDict, parDict, diagrams, newplot(), describe()
#------------------------------------------------------------------------------------------------------------------
   
# Create stateDict that later will be used to store final state and used for initialization in 'cont':
global stateDict; stateDict =  {}
stateDict = {variable.derivative.name:None for variable in model_description.modelVariables \
                                            if variable.derivative is not None}
stateDict.update(timeDiscreteStates) 

global stateDictInitial; stateDictInitial = {}
for key in stateDict.keys():
    if not key[-1] == ']':
         if key[-3:] == 'I.y':
            stateDictInitial[key] = key[:-10]+'I_start'
         elif key[-3:] == 'D.x':
            stateDictInitial[key] = key[:-10]+'D_start'
         else:
            stateDictInitial[key] = key+'_start'
    elif key[-3] == '[':
        stateDictInitial[key] = key[:-3]+'_start'+key[-3:]
    elif key[-4] == '[':
        stateDictInitial[key] = key[:-4]+'_start'+key[-4:]
    elif key[-5] == '[':
        stateDictInitial[key] = key[:-5]+'_start'+key[-5:] 
    else:
        print('The state vector has more than 1000 states')
        break

global stateDictInitialLoc; stateDictInitialLoc = {}
for value in stateDictInitial.values(): stateDictInitialLoc[value] = value

# Create parDict
global parDict; parDict = {}
parDict['V_start']    = 0.35          # L
parDict['VXv_start'] = 0.35*0.2       
parDict['VXd_start'] = 0.0            
parDict['VG_start'] = 0.35*18.0       
parDict['VGn_start'] = 0.35*2.4       
parDict['VL_start'] = 0.0             
parDict['VN_start'] = 0.0             

parDict['qG_max1'] = 0.2971
parDict['qG_max2'] = 0.0384
parDict['qGn_max1'] = 0.1238
parDict['qGn_max2'] = 0.0218
parDict['mu_d_max'] = 0.1302

parDict['k_lysis'] = 0.0

parDict['alpha'] = -1.0
parDict['beta'] = 0.01

parDict['G_in']  =  15.0          # mM
parDict['Gn_in']  =  4.0          # mM
parDict['t0'] =   0.0             # h
parDict['F0'] =   0.0             # L/h
parDict['t1'] =  60.0             # h
parDict['F1'] =   0.00            # L/h
parDict['t2'] = 992.0             # h
parDict['F2'] =   0.00            # L/h
parDict['t3'] = 993.0             # h
parDict['F3'] =   0.0040          # L/h
parDict['t4'] = 994.0             # h
parDict['F4'] =   0.0080          # L/h
parDict['t5'] = 995.0             # h
parDict['F5'] =   0.012           # L/h
parDict['t6'] = 996.0             # h
parDict['F6'] =   0.012           # L/h

global parLocation; parLocation = {}
parLocation['V_start'] = 'bioreactor.V_start'
parLocation['VXv_start'] = 'bioreactor.m_start[1]'
parLocation['VXd_start'] = 'bioreactor.m_start[2]'
parLocation['VG_start'] = 'bioreactor.m_start[3]'
parLocation['VGn_start'] = 'bioreactor.m_start[4]'
parLocation['VL_start'] = 'bioreactor.m_start[5]'
parLocation['VN_start'] = 'bioreactor.m_start[6]'

parLocation['qG_max1'] = 'bioreactor.culture.qG_max1'
parLocation['qG_max2'] = 'bioreactor.culture.qG_max2'
parLocation['qGn_max1'] = 'bioreactor.culture.qGn_max1'
parLocation['qGn_max2'] = 'bioreactor.culture.qGn_max2'
parLocation['mu_d_max'] = 'bioreactor.culture.mu_d_max'

parLocation['alpha'] = 'bioreactor.culture.alpha'
parLocation['beta'] = 'bioreactor.culture.beta'

parLocation['k_lysis'] = 'bioreactor.broth_decay.k_lysis'

parLocation['G_in'] = 'feedtank.c_in[3]'
parLocation['Gn_in'] = 'feedtank.c_in[4]'

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
global key_variables; key_variables = []
parLocation['mu'] = 'bioreactor.culture.mu'; key_variables.append(parLocation['mu'])
parLocation['mu_d'] = 'bioreactor.culture.mu_d'; key_variables.append(parLocation['mu_d'])

# Parameter value check - especially for hysteresis to avoid runtime error
global parCheck; parCheck = []
parCheck.append("parDict['V_start'] > 0")
parCheck.append("parDict['VXv_start'] >= 0")
parCheck.append("parDict['VG_start'] >= 0")
parCheck.append("parDict['VGn_start'] >= 0")
parCheck.append("parDict['VL_start'] >= 0")
parCheck.append("parDict['VN_start'] >= 0")
parCheck.append("parDict['t0'] <= parDict['t1']")
parCheck.append("parDict['t1'] <= parDict['t2']")
parCheck.append("parDict['t2'] <= parDict['t3']")
parCheck.append("parDict['t3'] <= parDict['t4']")
parCheck.append("parDict['t4'] <= parDict['t5']")
parCheck.append("parDict['t5'] <= parDict['t6']")

# Create list of diagrams to be plotted by simu()
global diagrams
diagrams = []

def newplot(title='Fedbatch cultivation',  plotType='TimeSeries'):
   """ Standard plot window,
        title = '' """
   
   # Reset pens
   setLines()

   # Globals
   global ax11, ax12, ax21, ax22, ax31, ax32, ax41, ax42, ax51, ax52  
   global ax13, ax23, ax33, ax43
    
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
      ax31.set_ylabel('Viable cell conc [1E6/mL]')

      ax32.grid()
      ax32.set_ylabel('Dead cell conc [1E6/mL]')

      ax41.grid()
      ax41.set_ylabel('Feed rate [L/h]')
      ax41.set_xlabel('Time [h]')

      ax42.grid()
      ax42.set_ylabel('Volume [L]')
      ax42.set_xlabel('Time [h]')

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(sim_res['time'],sim_res['bioreactor.c[3]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(sim_res['time'],sim_res['bioreactor.c[5]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(sim_res['time'],sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(sim_res['time'],sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(sim_res['time'],sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(sim_res['time'],sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")       
      diagrams.append("ax41.plot(sim_res['time'],sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax42.plot(sim_res['time'],sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

   if plotType == 'Textbook_1':
 
      plt.figure()
      ax11 = plt.subplot(5,2,1); ax12 = plt.subplot(5,2,2)
      ax21 = plt.subplot(5,2,3); ax22 = plt.subplot(5,2,4)
      ax31 = plt.subplot(5,2,5); ax32 = plt.subplot(5,2,6)
      ax41 = plt.subplot(5,2,7)
      ax51 = plt.subplot(5,2,9) 

      ax11.set_title(title)
      ax11.grid()
      ax11.set_ylabel('Glucose [mM]')

      ax12.grid()
      ax12.set_ylabel('Lactate [mM]')

      ax21.grid()
      ax21.set_ylabel('Glutamine [mM]')

      ax22.grid()
      ax22.set_ylabel('Ammonia [mM]')

      ax31.grid()
      ax31.set_ylabel('Viable cell [1E6/mL]')

      ax32.grid()
      ax32.set_ylabel('Dead cell conc [1E6/mL]')
      ax32.set_xlabel('Time [h]')

      ax41.grid()
      ax41.set_ylabel('Feed rate [L/h]')
      
      ax51.grid()
      ax51.set_ylabel('Volume [L]')      
      ax51.set_xlabel('Time [h]')


      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(sim_res['time'],sim_res['bioreactor.c[3]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(sim_res['time'],sim_res['bioreactor.c[5]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(sim_res['time'],sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(sim_res['time'],sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(sim_res['time'],sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(sim_res['time'],sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")       
      diagrams.append("ax41.plot(sim_res['time'],sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax51.plot(sim_res['time'],sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

   if plotType == 'Textbook_2':
 
      plt.figure()
      ax11 = plt.subplot(5,3,1); ax12 = plt.subplot(5,3,2); ax13 = plt.subplot(5,3,3)
      ax21 = plt.subplot(5,3,4); ax22 = plt.subplot(5,3,5); ax23 = plt.subplot(5,3,6)
      ax31 = plt.subplot(5,3,7); ax32 = plt.subplot(5,3,8); ax33 = plt.subplot(5,3,9)
      ax41 = plt.subplot(5,3,10)
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
      
      ax51.grid()
      ax51.set_ylabel('Volume [L]')      
      ax51.set_xlabel('Time [h]')
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(sim_res['time'],sim_res['bioreactor.c[3]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(sim_res['time'],sim_res['bioreactor.c[5]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(sim_res['time'],sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(sim_res['time'],sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(sim_res['time'],sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(sim_res['time'],sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")       
      diagrams.append("ax41.plot(sim_res['time'],sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax51.plot(sim_res['time'],sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

      diagrams.append("ax13.set_title('- microscopic world')") 
      diagrams.append("ax13.plot(sim_res['time'],-(sim_res['bioreactor.culture.q[3]']+sim_res['bioreactor.culture.qG_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax13.plot(sim_res['time'],-sim_res['bioreactor.culture.q[3]'], color='b', linestyle=linetype)") 
      diagrams.append("ax23.plot(sim_res['time'],-(sim_res['bioreactor.culture.q[4]']+sim_res['bioreactor.culture.qGn_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax23.plot(sim_res['time'],-sim_res['bioreactor.culture.q[4]'], color='b', linestyle=linetype)") 
      diagrams.append("ax33.plot(sim_res['time'],sim_res['bioreactor.culture.q[1]'], color='b', linestyle=linetype)") 
      
      diagrams.append("ax11.set_ylim(0)")
      diagrams.append("ax13.set_ylim(0)")
      #diagrams.append("ax31.set_ylim(0)")
      diagrams.append("ax32.set_ylim(ax31.get_ylim())")
      #diagrams.append("plt.show()")

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
      diagrams.append("ax11.plot(sim_res['time'],sim_res['bioreactor.c[3]'], color='b', linestyle=linetype)")       
      diagrams.append("ax12.plot(sim_res['time'],sim_res['bioreactor.c[5]'], color='r', linestyle=linetype)")       
      diagrams.append("ax21.plot(sim_res['time'],sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax22.plot(sim_res['time'],sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax31.plot(sim_res['time'],sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax32.plot(sim_res['time'],sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")       
      diagrams.append("ax41.plot(sim_res['time'],sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)") 
      diagrams.append("ax42.plot(sim_res['time'],sim_res['bioreactor.c[7]'], color='g', linestyle=linetype)")       
      diagrams.append("ax51.plot(sim_res['time'],sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

      diagrams.append("ax13.set_title('- cell specific rates')") 
      diagrams.append("ax13.plot(sim_res['time'],-(sim_res['bioreactor.culture.q[3]']+sim_res['bioreactor.culture.qG_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax13.plot(sim_res['time'],-sim_res['bioreactor.culture.q[3]'], color='b', linestyle=linetype)") 
      diagrams.append("ax23.plot(sim_res['time'],-(sim_res['bioreactor.culture.q[4]']+sim_res['bioreactor.culture.qGn_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax23.plot(sim_res['time'],-sim_res['bioreactor.culture.q[4]'], color='b', linestyle=linetype)") 
      diagrams.append("ax33.plot(sim_res['time'],sim_res['bioreactor.culture.q[1]'], color='b', linestyle=linetype)") 
      diagrams.append("ax43.plot(sim_res['time'],sim_res['bioreactor.culture.q[7]'], color='g', linestyle=linetype)") 
      
      diagrams.append("ax11.set_ylim(0)")
      diagrams.append("ax13.set_ylim(0)")
      diagrams.append("ax32.set_ylim(ax31.get_ylim())")

      
def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""

   if name == 'culture':
      print('Reactor culture CHO-MAb - cell line HB-58 American Culture Collection ATCC') 

   elif name in ['broth', 'liquidphase', 'liquid-phase''media']:

      Xv  = model_get('liquidphase.Xv')[0]; 
      Xv_description = model_get_variable_description('liquidphase.Xv'); 
      Xv_mw = model_get('liquidphase.mw[1]')[0]
      
      Xd = model_get('liquidphase.Xd')[0]; 
      Xd_description = model_get_variable_description('liquidphase.Xd'); 
      Xd_mw = model_get('liquidphase.mw[2]')[0]
      
      G = model_get('liquidphase.G')[0]; 
      G_description = model_get_variable_description('liquidphase.G'); 
      G_mw = model_get('liquidphase.mw[3]')[0]
      
      Gn = model_get('liquidphase.Gn')[0]; 
      Gn_description = model_get_variable_description('liquidphase.Gn'); 
      Gn_mw = model_get('liquidphase.mw[4]')[0]
      
      L = model_get('liquidphase.L')[0]; 
      L_description = model_get_variable_description('liquidphase.L'); 
      L_mw = model_get('liquidphase.mw[5]')[0]
      
      N = model_get('liquidphase.N')[0]; 
      N_description = model_get_variable_description('liquidphase.N'); 
      N_mw = model_get('liquidphase.mw[6]')[0]
      
      Pr = model_get('liquidphase.Pr')[0]; 
      Pr_description = model_get_variable_description('liquidphase.Pr'); 
      Pr_mw = model_get('liquidphase.mw[7]')[0]

      print('Reactor broth substances included in the model')
      print()
      print(Xv_description, 'index = ', Xv, 'molecular weight = ', Xv_mw, 'Da')
      print(Xd_description, '  index = ', Xd, 'molecular weight = ', Xd_mw, 'Da')
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
FMU_explore = 'FMU-explore for FMPy version 1.0.0'
#------------------------------------------------------------------------------------------------------------------

# Define function par() for parameter update
def par(parDict=parDict, parCheck=parCheck, parLocation=parLocation, *x, **x_kwarg):
   """ Set parameter values if available in the predefined dictionaryt parDict. """
   x_kwarg.update(*x)
   x_temp = {}
   for key in x_kwarg.keys():
      if key in parDict.keys():
         x_temp.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an accessible parameter - check the spelling')
   parDict.update(x_temp)
   
   parErrors = [requirement for requirement in parCheck if not(eval(requirement))]
   if not parErrors == []:
      print('Error - the following requirements do not hold:')
      for index, item in enumerate(parErrors): print(item)

# Define function init() for initial values update
def init(parDict=parDict, *x, **x_kwarg):
   """ Set initial values and the name should contain string '_start' to be accepted.
       The function can handle general parameter string location names if entered as a dictionary. """
   x_kwarg.update(*x)
   x_init={}
   for key in x_kwarg.keys():
      if '_start' in key: 
         x_init.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an initial value, use par() instead - check the spelling')
   parDict.update(x_init)

# Define fuctions similar to pyfmi model.get(), model.get_variable_descirption(), model.get_variable_unit()
def model_get(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get() but returns just a value and not a list"""
   par_var = model_description.modelVariables
   for k in range(len(par_var)):
      if par_var[k].name == parLoc:
         try:
            if par_var[k].name in start_values.keys():
                  value = start_values[par_var[k].name]
            elif par_var[k].variability in ['constant', 'fixed']: 
               if par_var[k].type in ['Integer', 'Real']: 
                  value = float(par_var[k].start)      
               if par_var[k].type in ['String']: 
                  value = par_var[k].start                        
            elif par_var[k].variability == 'continuous':
               try:
                  timeSeries = sim_res[par_var[k].name]
                  value = timeSeries[-1]
               except (AttributeError, ValueError):
                  value = None
                  print('Variable not logged')
            else:
               value = None
         except NameError:
            print('Error: Information available after first simution')
            value = None          
   return value

def model_get_variable_description(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_description() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].description) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.description for x in par_var if parLoc in x.name]   
   return value[0]
   
def model_get_variable_unit(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_unit() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].unit) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.unit for x in par_var if parLoc in x.name]
   return value[0]
      
# Define function disp() for display of initial values and parameters
def disp(name='', decimals=3, mode='short'):
   """ Display intial values and parameters in the model that include "name" and is in parLocation list.
       Note, it does not take the value from the dictionary par but from the model. """
   
   def dict_reverser(d):
      seen = set()
      return {v: k for k, v in d.items() if v not in seen or seen.add(v)}
   
   if mode in ['short']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:
               print(dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
            else:
               print(dict_reverser(parLocation)[Location] , ':', model_get(Location))               
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parName,':', np.round(model_get(parLocation[parName]),decimals))
               else: 
                  print(parName,':', model_get(parLocation[parName])[0])

   if mode in ['long','location']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:       
               print(Location,':', dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parLocation[parName], ':', dict_reverser(parLocation)[Location], ':', parName,':', 
                     np.round(model_get(parLocation[parName]),decimals))

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

# Define simulation
def simu(simulationTime=simulationTime, mode='Initial', options=opts_std, diagrams=diagrams):
   """Model loaded and given intial values and parameter before, and plot window also setup before."""   
   
   # Global variables
   global sim_res, prevFinalTime, stateDict, stateDictInitial, stateDictInitialLoc, start_values
   
   # Simulation flag
   simulationDone = False
   
   # Internal help function to extract variables to be stored
   def extract_variables(diagrams):
       output = []
       variables = [v for v in model_description.modelVariables if v.causality == 'local']
       for j in range(len(diagrams)):
           for k in range(len(variables)):
               if variables[k].name in diagrams[j]:
                   output.append(variables[k].name)
       return output

   # Run simulation
   if mode in ['Initial', 'initial', 'init']: 
      
      start_values = {parLocation[k]:parDict[k] for k in parDict.keys()}
      
      # Simulate
      sim_res = simulate_fmu(
         filename = fmu_model,
         validate = False,
         start_time = 0,
         stop_time = simulationTime,
         output_interval = simulationTime/options['ncp'],
         record_events = True,
         start_values = start_values,
         fmi_call_logger = None,
         output = list(set(extract_variables(diagrams) + list(stateDict.keys()) + key_variables))
      )
      
      simulationDone = True
      
   elif mode in ['Continued', 'continued', 'cont']:
      
      if prevFinalTime == 0: 
         print("Error: Simulation is first done with default mode = init'")
         
      else:         
         # Update parDictMod and create parLocationMod
         parDictRed = parDict.copy()
         parLocationRed = parLocation.copy()
         for key in parDict.keys():
            if parLocation[key] in stateDictInitial.values(): 
               del parDictRed[key]  
               del parLocationRed[key]
         parLocationMod = dict(list(parLocationRed.items()) + list(stateDictInitialLoc.items()))
   
         # Create parDictMod and parLocationMod
         parDictMod = dict(list(parDictRed.items()) + 
            [(stateDictInitial[key], stateDict[key]) for key in stateDict.keys()])      

         start_values = {parLocationMod[k]:parDictMod[k] for k in parDictMod.keys()}
  
         # Simulate
         sim_res = simulate_fmu(
            filename = fmu_model,
            validate = False,
            start_time = prevFinalTime,
            stop_time = prevFinalTime + simulationTime,
            output_interval = simulationTime/options['ncp'],
            record_events = True,
            start_values = start_values,
            fmi_call_logger = None,
            output = list(set(extract_variables(diagrams) + list(stateDict.keys()) + key_variables))
         )
      
         simulationDone = True
   else:
      
      print("Error: Simulation mode not correct")

   if simulationDone:
      
      # Plot diagrams from simulation
      linetype = next(linecycler)    
      for command in diagrams: eval(command)
   
      # Store final state values in stateDict:        
      for key in stateDict.keys(): stateDict[key] = model_get(key)  
         
      # Store time from where simulation will start next time
      prevFinalTime = sim_res['time'][-1]
      
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
    
#   variables = list(model.get_model_variables().keys())
   variables = [v.name for v in model_description.modelVariables]
        
   for i in range(len(variables)):
      component = model_component(variables[i])
      if (component not in component_list) \
      & (component not in ['','BPL', 'Customer', 'today[1]', 'today[2]', 'today[3]', 'temp_2', 'temp_3']):
         component_list.append(component)
      
   print(sorted(component_list, key=str.casefold))

# Describe MSL   
def describe_MSL(flag_vendor=flag_vendor):
   """List MSL version and components used"""
   print('MSL:', MSL_usage)
 
# Describe parameters and variables in the Modelica code
def describe_general(name, decimals):
  
   if name == 'time':
      description = 'Time'
      unit = 'h'
      print(description,'[',unit,']')
      
   elif name in parLocation.keys():
      description = model_get_variable_description(parLocation[name])
      value = model_get(parLocation[name])
      try:
         unit = model_get_variable_unit(parLocation[name])
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
      description = model_get_variable_description(name)
      value = model_get(name)
      try:
         unit = model_get_variable_unit(name)
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
   print('Model for bioreactor has been setup. Key commands:')
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
#   FMU_type = model.__class__.__name__
   constants = [v for v in model_description.modelVariables if v.causality == 'local']
   
   print()
   print('System information')
   print(' -OS:', platform.system())
   print(' -Python:', platform.python_version())
   try:
       scipy_ver = scipy.__version__
       print(' -Scipy:',scipy_ver)
   except NameError:
       print(' -Scipy: not installed in the notebook')
   print(' -FMPy:', version('fmpy'))
   print(' -FMU by:', read_model_description(fmu_model).generationTool)
   print(' -FMI:', read_model_description(fmu_model).fmiVersion)
   if model_description.modelExchange is None:
      print(' -Type: CS')
   else:
      print(' -Type: ME')
   print(' -Name:', read_model_description(fmu_model).modelName)
   print(' -Generated:', read_model_description(fmu_model).generationDateAndTime)
   print(' -MSL:', MSL_version)    
   print(' -Description:', BPL_version)   
   print(' -Interaction:', FMU_explore)
   
#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

BPL_info()