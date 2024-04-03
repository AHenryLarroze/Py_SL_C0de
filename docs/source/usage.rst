Usage
=====


INSTALLATION
------------

.. note:: The code is developped on windows for the moment. We want to export it on other OS but for the moment we focus on the packages. 

Download 

To install the library, run in command line : 

.. code-block::
    
    (.venv) $ conda create -n SL_C0de 
    (.venv) $ conda activate SL_C0de
    (.venv) $ conda install pip      
    (.venv) $ pip install --index-url https://test.pypi.org/simple/ --no-deps slcode 

To install the dependencies, download the requirements.txt on https://github.com/AHenryLarroze/Py_SL_C0de and run inside the file :

.. code-block::

    (.venv) $ pip install -r requirements

You also need the stripy package. This package has unstable deployment. Try to install this package with pip and conda. 

.. note:: You need to be in the activated environment to install the packages locally

For testing the module. You can download the script file in https://github.com/AHenryLarroze/Py_SL_C0de. You will also need cartopy to install it. 

.. code-block::

    (.venv) $ conda install -c conda-forge cartopy


minimal exemple
---------------

Here we propose to run the code over 26 kyr with a time resolution of 500 years. The spatial resolution is deduced from the maximum spherical harmonic degree (see section xx). Here 64 degree correspond to a spatial resolution of 170 km at the equator. 

.. code-block::

    import numpy as np

    stop=26 #stop the computation at 26 kyr
    step=0.5 # run the model at a time resolution of 0.5
    time_step=np.arange(start=stop,stop=-step,step=-step)
    maxdeg=64 #Define the maximum degrre of spherical harmonics

We preset the ice, sediment and topographic time grid.

.. code-block::

    from SL_C0de.SOLVER import Precomputation

    from scipy import io
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    # Loading Ice data

    data=io.loadmat(r'H:/SL_C0de/SL_C0de_data/ice6G122k') #load the file, can be modified to load an other file.
    ice_in=data['ice6g']# ice_in[temps, latitude, longitude] hauteur de la glace
    ice_time=data['ice_time']
    ice_time=ice_time.squeeze()
    ice_time=np.append(ice_time[0]+1,ice_time)
    ice_lon =data['ice_long'].squeeze()
    ice_lat =data['ice_lat'].squeeze()[::-1]
    # The ice is not in the good shape. We have to derivate it : 
    ice_0=ice_in.T[0,:,:].copy()
    ice_in=np.diff(ice_in.T,axis=0)
    ice_in=np.concatenate((np.expand_dims(ice_0,axis=0),ice_in),axis=0)

    ice_grid=dict(name='ICE_ICE6G', grid=ice_in, time=ice_time,lon=ice_lon,lat=ice_lat)

    # Loading sediment loading data

    sed_ncgrid = pd.read_csv(r'H:/SL_C0de/SL_C0de_data/AYS2_sed.csv') 

    lon=np.array(sed_ncgrid["POINT_X.N.19.11"])
    lat=np.array(sed_ncgrid["POINT_Y.N.19.11"])
    sed=np.array(sed_ncgrid[:]).T/4
    sed=sed[3:,:]
    sed=sed[::-1,:]*0 # add a zero sediment model
    sed_time_step=np.arange(26.25,-0.25,-0.25)

    area=(93,98,13,20)

    sed_grid=dict(name='SED', grid=sed, time=sed_time_step,lon=lon,lat=lat,frame=area)

    # Loading topographic data

    data=io.loadmat(r'H:/SL_C0de/SL_C0de_data/topo_SL') #load the file, can be modified to load an other file.
    topo_grid=data['topo_bed'][::2,1::2].squeeze()
    topo_lat=data['lat_topo'].squeeze()[::-2]
    topo_lon=data['lon_topo'].squeeze()[1::2]

    topo_grid=dict(name='topo_SL',grid=topo_grid,lon=topo_lon,lat=topo_lat)

For this example we will set no ice load, with the method zeros_time. topography have a general depth of 3000 m. We create a disk centered at the equator with a thickness of 100 m, a radius of 300 km and deposited after 500 years. 

.. code-block::

    # Precompute the full data to prepare the next resolution of the SLE

    stop=26 # define the number of time steps
    step=0.5
    time_step=np.arange(start=stop,stop=-step,step=-step) # in kyr to present
    maxdeg=512 # define the maximum degree for spherical harmonics.
    Output_way='H:/SL_C0de/Interpolated_grid/Test_512_no_sed'
    plot=True

    Precomputation(ice_grid,sed_grid,topo_grid,Output_way,stop=stop,step=step,maxdeg=maxdeg,plot=plot)


Now all entries are set up, we can resolve the sea level equation. We set as entry love numbers based on VM5 :cite:p:`peltier_2015` with a litosphere of 60 km, a upper mantle viscosity of :math:`10^{20.5}` Pa.s and a lower mantle viscosity of :math:`10^{22.699}` Pa.s. We use a forward modeling version of the SLE resolution. We set a convergence criterion of the SLE to :math:`10^{-10}`

.. code-block::

    from SL_C0de.SOLVER import SLE_solver

    Input_way='H:/SL_C0de/Interpolated_grid/Test_512_no_sed'
    ice_name='ICE_ICE6G'
    sediment_name='AYS2'
    topo_name='topo_SL'
    ocean_name='OCE'
    love_way='H:/SL_C0de/SL_C0de_data/variable_um_32'
    love_file='VM5a.l32.um21.5.lm22.699'
    conv_lim=10e-10
    topo_lim=3
    Output_way='H:/SL_C0de/Interpolated_grid/Test_512_no_sed'

    load=SLE_solver(Input_way,ice_name,sediment_name,topo_name,ocean_name,love_way,love_file,topo_lim,conv_lim,Output_way)

From here, we can post process these data to compute subsidence linked to the different loads. First you must save your result into a file.

.. code-block::

    import os
    import sys

    Output_way="outputs/"

    if not(os.path.exists(Output_way+love_file)):
        os.mkdir(Output_way+love_file)
    ocean_time_grid.timecoefftotimegrd()
    ocean_time_grid.save(Output_way+love_file)
    ice_time_grid_model.save(Output_way+love_file)
    topo_time_grid_model.save(Output_way+love_file)
    sediment_time_grid.save(Output_way+love_file)

Now we apply the post process 

.. code-block::

    from SL_C0de.SOLVER import Post_process

    Input_way='H:/SL_C0de/Interpolated_grid/Test_512_no_sed/model_output'
    sed_name='AYS2'
    ice_name='ICE_ICE6G'
    ocean_name='OCE'
    topo_name='topo_SL'
    love_way='H:/SL_C0de/SL_C0de_data/variable_um_32'
    Post_process(Input_way,sed_name,ice_name,ocean_name,topo_name,love_way)

Post process create new files in a folder called LOAD in the output folder. These data can be plotted using different functions.

.. code-block::

    import numpy as np
    Input_way='H:/SL_C0de/Interpolated_grid/Test_512/model_output/VM5a.l32.um20.25.lm22.699'
    plot=dict(plot=True,
            times=[0,5,10,15,20],
            frames=[(-179,179,-89,89)],
            frames_resolution=[1024],
            frames_min_max=np.array([[None,None,None,None],[None,None,None,None]]),
            contours_v=[[[i for i in range(-30,-10,10)]+[i for i in range(-10,10,1)]+[i for i in range(10,30+1,10)]+[100,200],],[[-3,-0.1,0,0.1,3,10,20]]]
            transects=[(21.63,96.01,11.1235,96.04)],
            point_density=[512],
            transect_min_max=[(-2,0.5,-0.7,1),(-3,2,-0.7,1)],
            points=np.array([[16.8,96],[15,96.25],[15.9609036,95.4758372],[16.0318760,94.8894641],[16.3469326,95.2544328]]),
            )

.. code-block::

    from SL_C0de.SOLVER import plot_model_result_map,plot_model_result_cross_section,plot_model_output_points
    plot_model_result_map(Input_way,plot)
    plot_model_result_cross_section(Input_way,plot)
    plot_model_output_points(Input_way,plot)

.. code-block::

    from SL_C0de.SOLVER import export_to_netcdf
    export_to_netcdf(Input_way,0,plot['frames_resolution'][1],plot['frames'][1],Input_way+'\LOAD','SEDIMENT','LOAD')



