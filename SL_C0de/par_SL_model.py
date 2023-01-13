import numpy as np
def f_V_lm(beta_R_l,beta_G_l,t_it,sdelLw,sdelLi,sdelLs):
    V_lm_R_w=np.dot(beta_R_l.transpose()[t_it-1,:t_it-1].transpose(),sdelLw[:t_it-1].transpose())
    V_lm_R_i=np.dot(beta_R_l.transpose()[t_it-1,:t_it-1].transpose(),sdelLi[:t_it-1].transpose())
    V_lm_R_s=np.dot(beta_R_l.transpose()[t_it-1,:t_it-1].transpose(),sdelLs[:t_it-1].transpose())
    V_lm_G_w=np.dot(beta_G_l.transpose()[t_it-1,:t_it-1].transpose(),sdelLw[:t_it-1].transpose())
    V_lm_G_i=np.dot(beta_G_l.transpose()[t_it-1,:t_it-1].transpose(),sdelLi[:t_it-1].transpose())
    V_lm_G_s=np.dot(beta_G_l.transpose()[t_it-1,:t_it-1].transpose(),sdelLs[:t_it-1].transpose())
    return V_lm_R_w,V_lm_R_i,V_lm_R_s,V_lm_G_w,V_lm_G_i,V_lm_G_s

def f_V_lm_T(beta_tide,t_it,sdelLa):
    return np.dot(beta_tide.transpose()[t_it-1,:t_it-1].transpose(), sdelLa[:t_it-1].transpose())