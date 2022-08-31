import numpy as np
import numpy.matlib as npmlib

def f_love_lm(h,n):
    return npmlib.repmat(h,1,n+1) # reproduce the love number coefficient for n+1 time

def f_calc_beta(t_it,maxdeg,beta_l,beta_konly_l,mode_found,k_amp,h_amp,spoles,time_step):
    #initialise the temporary variable for calculation use
    beta_R_l_int=beta_l.copy()
    beta_R_konly_l_int=beta_konly_l.copy()
    beta_G_l_int=beta_l.copy()
    beta_G_konly_l_int=beta_konly_l.copy()
    for n in range(1,t_it):
        # initialize the beta for each time step
        beta_R = np.zeros((maxdeg,)) 
        beta_G = np.zeros((maxdeg,)) 
        for lm in range(maxdeg):
            # num_mod is the number of love number needed for each degree-order
            num_mod = mode_found[lm][0] 
            # calculate the beta coefficient for each time step and degree-order
            beta_R[lm]=np.sum((h_amp[lm,:num_mod])/spoles[lm,:num_mod] * (1 - np.exp(- spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
            beta_G[lm]=np.sum((k_amp[lm,:num_mod])/spoles[lm,:num_mod] * (1 - np.exp(- spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
        # add to the beta a coefficient of value 0 !!!!
        beta_R_l_int[n-1,:]=np.concatenate((np.array([0]),beta_R))
        beta_G_l_int[n-1,:]=np.concatenate((np.array([0]),beta_G))
        # for rotation only needed for degree 2
        lm=1
        num_mod=mode_found[lm][0]
        # calculate the beta konly 
        beta_G_konly_l_int[n-1] = np.sum((k_amp[lm,:num_mod])/spoles[lm,:num_mod] * (1 - np.exp(- spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
        beta_R_konly_l_int[n-1] = np.sum((h_amp[lm,:num_mod])/spoles[lm,:num_mod] * (1 - np.exp(- spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
    #for each time step update the beta coeffcient
    return beta_G_l_int,beta_G_konly_l_int,beta_R_l_int,beta_R_konly_l_int

def f_calc_beta_tide(t_it,maxdeg,beta_tide,beta_konly_tide,mode_found,k_amp_tide,h_amp_tide,spoles,time_step):
        beta_G_tide_int=beta_tide.copy()
        beta_G_konly_tide_int=beta_konly_tide.copy()
        beta_R_tide_int=beta_tide.copy()
        beta_R_konly_tide_int=beta_konly_tide.copy()
        for n in range(1,t_it):
            #initialize the beta coefficient at each time step time step for computational use
            beta_G = np.zeros((maxdeg,))
            beta_R = np.zeros((maxdeg,))
            for lm in range(maxdeg):
                num_mod = mode_found[lm][0]
                beta_G[lm]= sum((k_amp_tide[lm,:num_mod])/spoles[lm,:num_mod]* (1 - np.exp(- spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
                beta_R[lm]= sum((h_amp_tide[lm,:num_mod])/spoles[lm,:num_mod]* (1 - np.exp(- spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
            beta_G_tide_int[n-1,:]=np.concatenate((np.array([0]),beta_G))
            beta_R_tide_int[n-1,:]=np.concatenate((np.array([0]),beta_G))
            lm=1
            num_mod=mode_found[lm][0]
            beta_G_konly_tide_int[n-1] = np.sum((k_amp_tide[lm,:num_mod])/spoles[lm,:num_mod] * (1 - np.exp(- spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
            beta_R_konly_tide_int[n-1] = np.sum((h_amp_tide[lm,:num_mod])/spoles[lm,:num_mod] * (1 - np.exp(- spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
        return beta_G_tide_int,beta_G_konly_tide_int,beta_R_tide_int,beta_R_konly_tide_int