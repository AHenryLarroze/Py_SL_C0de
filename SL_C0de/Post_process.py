from SL_model import SL_model
from spharm import legendre, GaussQuad
from time import time
import numpy as np
if __name__ == '__main__':
    model=SL_model(maxdeg=255,nb_workers=12,topo_it_max=10,k_max=10) # Generate the model base
    model.set_from_file('C:/Users/ahenry01/Desktop/Python_code/SL_C0de_output/vl.VM5a.l40.um21.lm22.699tresult_25.09.2022_16.47.13')