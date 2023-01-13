import numpy as np

def f_V_lm(beta_l,t_it,sdelL):
    if t_it==1:
        return 0
    else :
        # return np.dot(beta_l[t_it-1,:t_it-1].transpose(),sdelL.transpose())
        return np.dot(beta_l.transpose()[t_it-1,:t_it-1].transpose(),sdelL[:t_it-1].transpose())
        
def f_V_lm_T(beta_tide,t_it,sdelLa):
    if t_it==1:
        return 0
    else :
        # return np.dot(beta_tide[t_it-1,:t_it-1].transpose(),sdelLa.transpose())
        return np.dot(beta_tide.transpose()[t_it-1,:t_it-1].transpose(), sdelLa[:t_it-1].transpose())
