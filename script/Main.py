from SL_C0de.SL_model import SL_model
from SL_C0de.spharm import legendre, GaussQuad
from time import time
import numpy as np



if __name__ == '__main__':
    model=SL_model(maxdeg=128,nb_workers=160,topo_it_max=10,k_max=10,data_way='C:/Users/ahenry01/Desktop/Python_code/SL_C0de_data',love_way='VM5a',time_step=np.array([26,25,24,23,22,21,20.5,20,19.5,19,18.5,18,17.5,17,16.5,16,15.5,15,14.5,14,13.5,13,12.5,12,11.5,11,10.5,10,9.5,9,8.5,8,7.5,7,6.5,6,5.5,5,4.5,4,3.5,3,2.5,2,1.5,1,0.5,0])) # Generate the model base
    #predefine spherical harmonics parameters

    model.P_lm=legendre(model.maxdeg,model.pool) #Calculate the Legendre associated functions for the Gaussian grid
    model.x, model.w= GaussQuad(model.maxdeg) # calculate the Gaussian grid parameters

    model.create_GRID() #create the grid framework

    model.create_parameters() #Create all grid parameters


    model.ice.quick_load(model.grid,model.pool,model.ice_way) #set the ice model
    model.sed.zeros(model) # set the sedimentation model
    model.topo.load(model,model.topo_way) # set the topography model, this function need the sediment and ice model ! 
    #Disk_load=model.grid.disk(model,0,180,1,100)
    # we need to calculate the ocean function after setting the topography due to it's dependencies
    model.create_ocean(model.topo.topo_pres) # create the ocean function (see theoria (TABOO), spada et al., 2003) from the actual topography.

    #regarder si il n'y a pas une boucle while plutot avec un critère de convergence pour stopper la boucle.
    # initialization of the prev variable
    for model.topo_it in range(model.topo_it_max): # Il y a un pb sur le chargement de topo à corriger !
        model.reset() #reset modify the variable for the next itteration but it also create some variable. Therefor, i need to check better the code

        model.ice.deli_00_prev = 0 #modifier cette entrès quand cette variable seras dans l'objet
        
        #set topography
        model.topo.topo[0,:,:]=model.topo.topo_initial.copy()      
        model.topo.ice_correction(model)
        model.topo.update_topo_0()
        
        
        model.oc.evaluate_ocean(model.topo.topo_0).grdtocoeff(model)
        model.oc.update_oc_0()
        model.oc.save_prev()
        model.SL.saved=np.array([])
        model.oc.saved=np.array([])

        model.reset_save_data()
        
        for t_it in range(1,len(model.grid.time_step)):
            t1=time()
            # initial topo with the ice correction : 
            model.topo.modify(model.topo.topo[t_it,:,:])
            # grd correspond donc au topo_j défini dans le code de kendal et al.
            model.oc.evaluate_ocean(model.topo.grd).grdtocoeff(model)
            model.SL.TO.modify(model.topo.topo_0*(model.oc.grd-model.oc.oc_0_grd)).grdtocoeff(model)
            
            model.sed.calc_del_sed(t_it).grdtocoeff(model)
            model.ice.calc_del_ice(t_it).grdtocoeff(model)
            
            model.ice.sdeli_00=model.ice.coeff[0]-model.ice.deli_00_prev # peut être à passer dans l'initialisation ou le modify de ice object après leurs création
            
            #apply the loop solving the see level equation
            k = 0
            chi = model.epsilon * 2
            while (k < model.k_max) and (chi >= model.epsilon):
                model.SL.calc_sdelS(model,k)# calculate the small variation of sea level
                #Disk_load.modify(Disk_load.disk[t_it,:,:]-Disk_load.disk[t_it-1,:,:]).grdtocoeff(model)
                model.load.modify(t_it,model.ice.rho*model.ice.coeff + model.oc.rho*model.SL.delS.coeff+ model.sed.rho*model.sed.coeff)
                # calculate viscous contribution
                # beta contains the viscous love numbers for time t_it,
                # row index goes over the time increments, column
                # index goes over lm
                model.load.calc_viscuous_load(model,t_it,model.load.sdelL,model.love.beta_l)
                if model.rotation : 
                    model.load.calc_rotational_potential(model,t_it)
                model.load.calc_viscuous_load_T(model,t_it,model.load.sdelLa)
                # calculate sea level perturbation
                # add ice and sea level and multiply with love numbers
                model.SL.calc_sdelS_new(model)
                chi = np.abs((np.sum(np.abs(model.SL.sdelS_new)) - np.sum(np.abs(model.SL.delS.coeff))) / np.sum(np.abs(model.SL.delS.coeff)))
                model.SL.sdelS.modify(model.SL.sdelS_new.copy(),'coeff')
                k+=1
            #update the prev variable for the next loop
            model.save_data()
            model.save_prev()
            model.SL.ESL[t_it]=np.real(model.ice.coeff[0]/model.oc.area*model.ice.rho/model.oc.rho)

            #SL.ESL[t_it]=ice.deli_00_prev/oc.area*ice.rho/oc.rho
            model.SL.delSL.modify(model.SL.delSLcurl.grd+model.SL.delPhi_g)
            model.topo.topo[t_it,:,:]=-model.SL.delSL.grd+model.topo.topo_0#update cette ligne avec le nouvel objet spherical_top
            t2=time()
            print('time : ',model.grid.time_step[t_it],' in ',t2-t1)
        model.ice.topo_pres_ice_corrected = model.topo.topo_pres - model.ice.ice[-1,:,:] + model.ice.ice_corrected[-1,:,:]
    model.topo.topo_initial=model.ice.topo_pres_ice_corrected - (model.topo.topo[-1,:,:]-model.topo.topo[0,:,:]) 
    model.save_model(output_way=r'C:/Users/ahenry01/Desktop/Python_code/SL_C0de_output',subname='VM5a')
