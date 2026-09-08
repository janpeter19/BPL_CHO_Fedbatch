# setup data CHO_Fedbatch_fmpy 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-09-08 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt 
from fmpy import simulate_fmu
from fmpy import read_model_description

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
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
      print('Linux - run FMU pre-compiled OpenModelica') 
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
   opts_std = {'NCP': 500}
elif flag_type in ['ME', 'me']:
   opts_std = {'NCP': 500}
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
   MSL_usage = '4.1.0 - used components: RealInput, RealOutput, CombiTimeTable, Types' 
   MSL_version = '4.1.0'
   BPL_version = 'Bioprocess Library version 2.3.2' 
else:    
   print('There is no FMU for this platform')

#------------------------------------------------------------------------------------------------------------------
   
# Simulation time
simulationTime = 120.0
prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = {variable.derivative.name:None for variable in model_description.modelVariables \
                                            if variable.derivative is not None}
stateValue.update(timeDiscreteStates) 

stateValueInitial = {}
for key in stateValue.keys():
    if not key[-1] == ']':
         if key[-3:] == 'I.y':
            stateValueInitial[key] = key[:-10]+'I_start'
         elif key[-3:] == 'D.x':
            stateValueInitial[key] = key[:-10]+'D_start'
         else:
            stateValueInitial[key] = key+'_start'
    elif key[-3] == '[':
        stateValueInitial[key] = key[:-3]+'_start'+key[-3:]
    elif key[-4] == '[':
        stateValueInitial[key] = key[:-4]+'_start'+key[-4:]
    elif key[-5] == '[':
        stateValueInitial[key] = key[:-5]+'_start'+key[-5:] 
    else:
        print('The state vector has more than 1000 states')
        break

stateValueInitialLoc = {}
for value in stateValueInitial.values(): stateValueInitialLoc[value] = value

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture', 'bioreactor.broth_decay']

# Provide process diagram on disk
fmu_process_diagram ='BPL_CHO_Fedbatch_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, ax, lines
#------------------------------------------------------------------------------------------------------------------
   
# Create parValue
parValue = {}

# Create dictionary of parameters full name - to get values in the notebook
parLocation = {}

# Extra only for describe()
keyVariables = []
parLocation['mu'] = 'bioreactor.culture.mu'; keyVariables.append(parLocation['mu'])
parLocation['mu_d'] = 'bioreactor.culture.mu_d'; keyVariables.append(parLocation['mu_d'])
parLocation['bioreactor.c[3]'] = 'bioreactor.c[3]'; keyVariables.append(parLocation['bioreactor.c[3]'])

# Parameter value check - especially for hysteresis to avoid runtime error
parCheck = []

# Create list of diagrams to be plotted by simu()
diagrams = []

# Create an empty list axes to be defined in newplot() and plotted by simu() or show()
ax = []

# Create list of pens for the diagrams
lines = ['-','--',':','-.']

