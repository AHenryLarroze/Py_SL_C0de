import pyshtools as pysh
import math
import numpy as np
import numpy.matlib as npmlib

def f_legendre(N,x):
    P_lm_shtools=pysh.legendre.legendre(N,x)/2
    P_lm_shtools[0,0]=P_lm_shtools[0,0]*math.sqrt(2)
    return P_lm_shtools

def f_grdtocoeff(l,P_lm,w,F_ym):
    return np.sum(P_lm[:l+1,:].transpose()*npmlib.repmat(np.expand_dims(w,axis=0).transpose(),1,l+1)*F_ym[:,:l+1],axis=0)

def f_coefftogrd(l,P_lm,coeff,N):
    return  P_lm[:l+1,:].transpose()*npmlib.repmat(get_coeffs(coeff,l),N,1)

def get_coeffs(a_lm,n):
    if n == 0 :
        a_n = a_lm[0]
    else :
        gauss_sum = int(n*(n+1)/2)
        # position of current order in a_lm vector
        a_n = a_lm[gauss_sum:gauss_sum+n+1]
    return a_n

