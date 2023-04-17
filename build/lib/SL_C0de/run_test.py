from SL_model import SL_model
model_p=SL_model(maxdeg=16,nb_workers=32)
model_p.run()
model_p.init_resolution_plot(256)
model_p.plot_RSL_3D(0.5)
