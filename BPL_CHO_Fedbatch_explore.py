# setup application functions BPL_CHO_Fedbatch, dependent on previous import of functions from fmu_explore 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-08-21 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Specific application functions: newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

def newplot(title='Fedbatch cultivation',  plotType='TimeSeries'):
   """ Standard plot window,
        title = '' """
   
   # Reset pens
   resetPen()
   
   # Plot diagram 
   if plotType == 'TimeSeries':

      ax11 = plt.subplot(4,2,1); ax12 = plt.subplot(4,2,2)
      ax21 = plt.subplot(4,2,3); ax22 = plt.subplot(4,2,4)
      ax31 = plt.subplot(4,2,5); ax32 = plt.subplot(4,2,6)
      ax41 = plt.subplot(4,2,7); ax42 = plt.subplot(4,2,8)
      
      ax.clear()
      ax.append(ax11)
      ax.append(ax12)
      ax.append(ax21)
      ax.append(ax22) 
      ax.append(ax31)
      ax.append(ax32)
      ax.append(ax41)
      ax.append(ax42)           
        
      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('Glucose conc [mM]')

      ax[1].grid()
      ax[1].set_ylabel('Lactate conc [mM]')

      ax[2].grid()
      ax[2].set_ylabel('Glutamine conc [mM]')

      ax[3].grid()
      ax[3].set_ylabel('Ammonia conc [mM]')

      ax[4].grid()
      ax[4].set_ylabel('Viable cells conc [1E6/mL]')

      ax[5].grid()
      ax[5].set_ylabel('Dead cells conc [1E6/mL]')

      ax[6].grid()
      ax[6].set_ylabel('Feed rate [L/h]')
      ax[6].set_xlabel('Time [h]')

      ax[7].grid()
      ax[7].set_ylabel('Volume [L]')
      ax[7].set_xlabel('Time [h]')

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[1].plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[3].plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[4].plot(t,sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[5].plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")        
      diagrams.append("ax[6].plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax[7].plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 
   
   if plotType == 'TimeSeries1':

      plt.figure()
      ax11 = plt.subplot(4,2,1); ax[1] = plt.subplot(4,2,2)
      ax21 = plt.subplot(4,2,3); ax22 = plt.subplot(4,2,4)
      ax31 = plt.subplot(4,2,5); ax32 = plt.subplot(4,2,6)
      ax41 = plt.subplot(4,2,7); ax42 = plt.subplot(4,2,8)  

      ax.clear()
      ax.append(ax11)
      ax.append(ax12)
      ax.append(ax21)
      ax.append(ax22) 
      ax.append(ax31)
      ax.append(ax32)
      ax.append(ax41)
      ax.append(ax42)       

      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('Glucose conc [mM]')

      ax[1].grid()
      ax[1].set_ylabel('Lactate conc [mM]')

      ax[2].grid()
      ax[2].set_ylabel('Glutamine conc [mM]')

      ax[3].grid()
      ax[3].set_ylabel('Ammonia conc [mM]')

      ax[4].grid()
      ax[4].set_ylabel('Viable cells conc [1E6/mL]')

      ax[5].grid()
      ax[5].set_ylabel('Dead and lysed cells [1E6/mL]')

      ax[6].grid()
      ax[6].set_ylabel('Feed rate [L/h]')
      ax[6].set_xlabel('Time [h]')

      ax[7].grid()
      ax[7].set_ylabel('Volume [L]')
      ax[7].set_xlabel('Time [h]')

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[1].plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[3].plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[4].plot(t,sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[5].plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")  
      diagrams.append("ax[5].plot(t,sim_res['bioreactor.c[3]'], color='k', linestyle=linetype)")       
      diagrams.append("ax[6].plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax[7].plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

    
   # Plot diagram 
   if plotType == 'TimeSeries2':
   
      plt.figure()
      ax11 = plt.subplot(4,2,1); ax12 = plt.subplot(4,2,2)
      ax21 = plt.subplot(4,2,3); ax22 = plt.subplot(4,2,4)
      ax31 = plt.subplot(4,2,5); ax32 = plt.subplot(4,2,6)
      ax41 = plt.subplot(4,2,7); ax42 = plt.subplot(4,2,8) 

      ax.clear()
      ax.append(ax11)
      ax.append(ax12)
      ax.append(ax21)
      ax.append(ax22) 
      ax.append(ax31)
      ax.append(ax32)
      ax.append(ax41)
      ax.append(ax42)        

      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('Glucose conc [mM]')

      ax[1].grid()
      ax[1].set_ylabel('Lactate conc [mM]')

      ax[2].grid()
      ax[2].set_ylabel('Glutamine conc [mM]')

      ax[3].grid()
      ax[3].set_ylabel('Ammonia conc [mM]')

      ax[4].grid()
      ax[4].set_ylabel('Viable cells [1E6]')

      ax[5].grid()
      ax[5].set_ylabel('Dead and lysed cells [1E6/mL]')

      ax[6].grid()
      ax[6].set_ylabel('Feed rate [L/h]')
      ax[6].set_xlabel('Time [h]')

      ax[7].grid()
      ax[7].set_ylabel('Volume [L]')
      ax[7].set_xlabel('Time [h]')

      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[1].plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[2].plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[3].plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[4].plot(t,sim_res['bioreactor.m[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[5].plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")  
      diagrams.append("ax[5].plot(t,sim_res['bioreactor.c[3]'], color='k', linestyle=linetype)")       
      diagrams.append("ax[6].plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)")       
      diagrams.append("ax[7].plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 
      
   if plotType == 'Textbook_3':
 
      plt.figure()
      ax11 = plt.subplot(5,3,1); ax12 = plt.subplot(5,3,2); ax13 = plt.subplot(5,3,3)
      ax21 = plt.subplot(5,3,4); ax22 = plt.subplot(5,3,5); ax23 = plt.subplot(5,3,6)
      ax31 = plt.subplot(5,3,7); ax32 = plt.subplot(5,3,8); ax33 = plt.subplot(5,3,9)
      ax41 = plt.subplot(5,3,10); ax42 = plt.subplot(5,3,11); ax43 = plt.subplot(5,3,12)
      ax51 = plt.subplot(5,3,13) 

      ax.clear()
      ax.append(ax11) # 0  
      ax.append(ax12) # 1 
      ax.append(ax13) # 2
      ax.append(ax21) # 3
      ax.append(ax22) # 4
      ax.append(ax23) # 5
      ax.append(ax31) # 6
      ax.append(ax32) # 7      
      ax.append(ax33) # 8
      ax.append(ax41) # 9            
      ax.append(ax42) #10  
      ax.append(ax43) #11  
      ax.append(ax51) #12  

      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('Glucose [mM]')

      ax[1].grid()
      ax[1].set_ylabel('Lactate [mM]')
      
      ax[2].grid()
      ax[2].set_ylabel('qG []')

      ax[3].grid()
      ax[3].set_ylabel('Glutamine [mM]')

      ax[4].grid()
      ax[4].set_ylabel('Ammonia [mM]')
      
      ax[5].grid()
      ax[5].set_ylabel('qGn []')

      ax[6].grid()
      ax[6].set_ylabel('Viable cell [1E6/mL]')

      ax[7].grid()
      ax[7].set_ylabel('Dead cell conc [1E6/mL]')
      ax[7].set_xlabel('Time [h]')
      
      ax[8].grid()
      ax[8].set_ylabel('mu [1/h]')
      ax[8].set_xlabel('Time [h]')

      ax[9].grid()
      ax[9].set_ylabel('Feed rate [L/h]')

      ax[10].grid()
      ax[10].set_ylabel('mAb []')

      ax[11].grid()
      ax[11].set_ylabel('qP []')

      ax[12].grid()
      ax[12].set_ylabel('Volume [L]')      
      ax[12].set_xlabel('Time [h]')
      
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['bioreactor.c[4]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[1].plot(t,sim_res['bioreactor.c[6]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[3].plot(t,sim_res['bioreactor.c[5]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[4].plot(t,sim_res['bioreactor.c[7]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[6].plot(t,sim_res['bioreactor.c[1]'], color='b', linestyle=linetype)")       
      diagrams.append("ax[7].plot(t,sim_res['bioreactor.c[2]'], color='r', linestyle=linetype)")       
      diagrams.append("ax[9].plot(t,sim_res['bioreactor.inlet[1].F'], color='b', linestyle=linetype)") 
      diagrams.append("ax[10].plot(t,sim_res['bioreactor.c[8]'], color='g', linestyle=linetype)")       
      diagrams.append("ax[12].plot(t,sim_res['bioreactor.V'], color='b', linestyle=linetype)") 

      diagrams.append("ax[2].set_title('- cell specific rates')") 
      diagrams.append("ax[2].plot(t,-(sim_res['bioreactor.culture.q[4]']+sim_res['bioreactor.culture.qG_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax[2].plot(t,-sim_res['bioreactor.culture.q[4]'], color='b', linestyle=linetype)") 
      diagrams.append("ax[5].plot(t,-(sim_res['bioreactor.culture.q[5]']+sim_res['bioreactor.culture.qGn_over']), color='r', linestyle=linetype)") 
      diagrams.append("ax[5].plot(t,-sim_res['bioreactor.culture.q[5]'], color='b', linestyle=linetype)") 
      diagrams.append("ax[8].plot(t,sim_res['bioreactor.culture.q[1]'], color='b', linestyle=linetype)") 
      diagrams.append("ax[11].plot(t,sim_res['bioreactor.culture.q[8]'], color='g', linestyle=linetype)") 
      
      diagrams.append("ax[0].set_ylim(0)")
      diagrams.append("ax[2].set_ylim(0)")
      diagrams.append("ax[7].set_ylim(ax[6].get_ylim())")      

      
def describe(name, model=model, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""

   if name == 'culture':
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
#  Startup
#------------------------------------------------------------------------------------------------------------------

FMU_explore_info()