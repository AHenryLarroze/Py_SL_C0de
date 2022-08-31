from SL_model import SL_model
model_p=SL_model(maxdeg=512,nb_workers=12,topo_it_max=10,k_max=10)
model_p.run()
model_p.init_resolution_plot(256)
model_p.calc_RSL()

model_p.calc_delR(model_p.SL.delS.saved*model_p.oc.rho,model_p.SL.sdelS.saved*model_p.oc.rho)


model_p.plot_3D_at_time(model_p.delR_e+model_p.delR_v,0.5)

from plot import write_saved_to_csv
saved=model_p.load.delRv['ice'].saved+model_p.load.delRv['water'].saved+model_p.load.delRe['ice'].saved+model_p.load.delRe['water'].saved+model_p.load.delRv_T['ice'].saved+model_p.load.delRv_T['water'].saved+model_p.load.delRe_T['ice'].saved+model_p.load.delRe_T['water'].saved
write_saved_to_csv(model_p,saved,'','del_R_tot')

## Sauvegarde des données
Data_dict={'g' : model_p.g,'M_e' : model_p.M_e,'omega' : model_p.omega,'k_hydro' : model_p.k_hydro, 'G' : model_p.G,'k_f' : model_p.k_f, 'a' : model_p.a, 'CminA' : model_p.CminA, 'C' : model_p.C, 'maxdeg' : model_p.maxdeg, 'rotation' : model_p.rotation, 'k_max' : model_p.k_max, 'epsilon' : model_p.epsilon, 'topo_it_max': model_p.topo_it_max, 'nb_workers' : model_p.nb_workers, 'time_step' : model_p.grid.time_step, 'sed_way' : 'Irrawady_sedimentationgrid', 'rho_sed' : model_p.sed.rho, 'ice_way' : 'ice6g_data', 'ice_corrected' : model_p.ice.ice_corrected, 'rho_ice' : model_p.ice.rho, 'delS' : model_p.SL.delS.saved, 'sdelS' : model_p.SL.sdelS.saved, 'rho_sea' : model_p.oc.rho,  'love_way' : 'VM5.l40.um316.lm5.1024/', 'topo_way' : 'topo_SL', 'topo' : model_p.topo.topo}

from datetime import datetime

# datetime object containing current date and time
now = datetime.now()

# dd/mm/YY H:M:S
dt_string = now.strftime("%d.%m.%Y_%H.%M.%S")

print(dt_string)

import pickle
with open('result_'+dt_string+'.pkl', 'wb') as f:
    pickle.dump(Data_dict, f)

from SL_model import SL_model
model_p=SL_model(maxdeg=512,nb_workers=12,topo_it_max=10,k_max=10)

model_p.set_from_file('result_25.07.2022_23.33.22')

import numpy as np

delL_sed=np.zeros((model_p.time_step_number,int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2)))+0j
sdelL_sed=np.zeros((model_p.time_step_number,int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2)))+0j
for t_it in range(1,model_p.time_step_number):
    print(model_p.grid.time_step[t_it])
    delL_sed[t_it,:]=model_p.sed.calc_del_sed(t_it).grdtocoeff(model_p).coeff.copy()
    if t_it==1:
        delL_sed_prev=delL_sed[t_it,:].copy()
    else :
        sdelL_sed[t_it-1,:]=delL_sed[t_it]-delL_sed_prev
        delL_sed_prev=delL_sed[t_it,:].copy()

model_p.calc_delR(delL_sed*model_p.sed.rho,sdelL_sed*model_p.sed.rho)

model_p.init_resolution_plot(512)
model_p.plot_map(model_p.delR_e+model_p.delR_v,0.5)
