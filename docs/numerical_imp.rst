Numerical implementation
========================

.. image:: /image/workflow.*

.. _iteration_desc:

The numerical implementation is using the derived impletation of :cite:p:`dalca_2013` from :cite:p:`kendall_2005`. This implementation is using two iteration counter, :math:`i` and :math:`j`. :math:`j` is associated to the time iteration and at each iteration the implementation iterate a convergence to find the best solution of the SLE where the counter is :math:`i`. The :math:`i` counter will be called the inner iteration. From :math:`\delta S^{i=1}_j` to :math:`\delta S^{i=\infty}_j` to will calculate the best solution of SLE. The numerical implementation require a third counter, :math:`k` called outer iteration, where the loop over the whole considered time is covered to improve first guess of initial topography, :math:`T_0^{k=1}` until the convergence (:math:`T_0^{k=\infty}`).

The SLE is rewrite to : 

.. math::
    \Delta \mathcal{SL}_j = \int_{-\infty}^{t_j} \iint_\Omega \Delta L (\Theta ^\prime, \Psi ^\prime,t ^\prime) \cdot [\frac{\Phi(\gamma, t_j - t^\prime)}{g} - \Gamma(\gamma,t_j-t^\prime)] \mathrm{d} \Omega ^\prime \mathrm{d} t^\prime - \Delta H_j - \Delta I_j 

Sea Level equation resolution
+++++++++++++++++++++++++++++

.. _SLE_res:

Sea level equation Resolution implementation
--------------------------------------------

We introduce the variation of ocean thickness : 

.. math::
    \delta S^{i,k}_j=-\Delta S^{i=\infty,k}_{j-1} + \Delta \mathcal{SL}^{i-1,k}_j C^{k-1}_j + \frac{\Delta \Phi (t_j)^{i-1,k}}{g} C^{k-1}_j - T_0^{k-1}[C_j^{k-1} - C_0^{k-1}]
    
If the last step of the outer iteration was completed, :math:`k-1`. The topography is updated :

.. math::
    T_j^{k-1} = T_p  + \Delta SL_p^{i=\infty,k-1} - \Delta SL_j^{i=\infty,k-1}

Where :math:`T_p` is the present day topography, we substract to the present day topography the earth movement. The resulting ocean function deduced from the reformulation of the topography is :

.. _oc_func:

.. math::
    C_j^{k-1} = 
    \begin{cases}
        1 & \text{if } T_j^{k-1} <0 \\
        0 & \text{if } T_j^{k-1} \geq 0
    \end{cases}

The sea level can be then estimated using a new estimation of the sea level change for the kth iteration.

.. math::
    \Delta \mathcal{SL}^{i-1,k}_j = \Delta \mathcal{G}^{i-1,k}_j - (\Delta R^{i-1,k}_j+\Delta H_j + \Delta I_j^{k-1})

for the first inner iteration :math:`i=1`, the initial variation of ocean thickness is predefine. The spatially invariant component is resulting from the variable ocean surface.

.. math::
    \frac{\Delta \Phi ^{i-1,k}_j}{g} = \frac{-1}{\mathcal{A}^{k-1}_j} \frac{\rho_I}{\rho_w}\iint_\Omega\Delta I_j^{k-1} d\Omega - \frac{-1}{\mathcal{A}^{k-1}_j} \iint_\Omega \Delta \mathcal{SL}_j^{i-1,k} C_j^{k-1} d\Omega + \iint_\Omega T_0^{k-1}[C_j^{k-1} - C_0^{k-1}] d\Omega

With 

.. math:: 
    \mathcal{A}^{k-1}_j=\iint_\Omega C_j^{k-1} d\Omega

Resolution of SLE including the deconvolution
------------------------------------------

.. _spec_sol:

The implementation in iteration result in a modification of the geoid and ground displacement :

.. math:: 
    \Delta \chi_j= \sum_{lm} T_l \sum_{n=0}^{j-1} \delta M_{lm} (t_n) Y_{lm}(\Theta,\Psi) \cdot I G F\left(\gamma, t_j-t_n \right)

This applied to the SLE equation, by linearity of the IGFs: 

.. math::
    [\Delta \mathcal{SL}_{lm,j}]^{i-1,k}=T_l E_l \Delta M_{lm,j}^{k,i} + T_l \sum_{n=0}^{j-1} \beta (l,t_n,t_j)\delta M_{lm,n}^{k,i} +\frac{1}{g}E^T_l([\Delta \Lambda_{lm,j-1}]^{i=\infty,k} + [\delta \Lambda_{lm,j}]^{i-1,k})+ \frac{1}{g} \sum^{j-1}_{n=0} \beta^T(l,t_n,t_j)[\delta \Lambda_{lm,j}]^{i=\infty,k} - \Delta H_{lm,n}-[\Delta I_{lm,n}]^{k-1}

where :math:`E_l = 1 + k_l^E - h_l^E`, :math:`\beta(l,t_n,t_j)=k_l^V(t_j-t_n)-h_l^V(t_j-t_n)` and :math:`T_l = \frac{4\pi a^2}{2l+1}`

In SL_C0de, :math:`T_l \sum_{n=0}^{j-1} \beta (l,t_n,t_j)\delta M_{lm,n}^{k,i}` is resolved in matrix produce. This result in a strong allocation of RAM as the :math:`\beta (l,t_n,t_j)` are stored in a matrix of size (time,time,(maximum degree + 1)(maximum degree +2)/2). The resulting time gain is very important. 

The conservation formula become :

.. math::
    \frac{\Delta \Phi_j}{g} = \frac{1}{C_{00,j}}(-\frac{\rho_i}{rho_w}\Delta I_{00,j}-RO_{00,j}+TO_{00,0})

Where :math:`RO_j = \Delta \mathcal{SL}_j C_j` and :math:`TO_j=T_0[C_j-C_0]`.

Convergence parameter
---------------------

Inner convergence on SLE
^^^^^^^^^^^^^^^^^^^^^^^^

.. _conv:

We define a convergence criterion :

.. math::
    \xi^{i,k}_j=|\frac{\sum_{l,m}|[\delta S_{lm}(t_j)]^{i,k}|-\sum_{l,m}|[\delta S_{lm}(t_j)]^{i-1,k}|}{\sum_{l,m}|[\delta S_{lm}(t_j)]^{i-1,k}|}|

Convergence for the SLE is limited by the convergence criterion : :math:`\xi_j^{i,k}`. We suppose that :math:`\xi_j^{i,k} < \epsilon_1` when 

.. math::
    [\delta S_{lm}(t_j)]^{i,k}=[\delta S_{lm}(t_j)]^{i=\infty,k}

Outer convergence criterion
^^^^^^^^^^^^^^^^^^^^^^^^^^^



Grounded ice correction
-----------------------

.. _Ice_corr:

The marine grounded ice is dependent of RSL variations. The ice is grounded if it satisfies : 

.. math::
    I_j > (SL_j + I_j)\frac{\rho_w}{\rho_I}

At each topographic iteration (:math:`k`) we update the grounded ice.

.. math::
    I_j^k = 
    \begin{cases}
        Ice\;Height & SL_j^{k-1} + Ice\;Height < 0 \\
        Ice\;Height & SL_j^{k-1} + Ice\;Height > 0 \\
         & and\;Ice\;Height > SL_j^{k-1},\frac{\rho_w}{\rho_I-\rho_w} \\
        0 & elsewhere
    \end{cases}

Computation of ground and geoid subsidence from different load source
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. _G_R_comp:

A functionality developed in :math:`SL_{C0de}` is the computation of the different component of the SLE separately, by type of Load and by viscous or elastic component. The development of this functionality was motivated by the necessity of exploring the different source of the RSL variation in a more and more complex modelization.

Elastic components of SLE :
---------------------------

We define 4 elastique component in the SLE, the ground displacement :math:`\Delta R^E_{lm}`, the geoid displacement :math:`\Delta G^E_{lm}`, the rotational ground displacement :math:`\Delta R^{T,E}_{lm}` and the rotational geoid displacement :math:`\Delta G^{T,E}_{lm}`.


.. math::
    \Delta R^E_{lm,j} = T_l h_l^E \Delta M_{lm,j}^{k,i}

.. math::
    \Delta G^E_{lm,j} = T_l (1+k_l^E) \Delta M_{lm,j}^{k,i}

.. math::
    \Delta R^{T,E}_{lm,j}=\frac{1}{g}h^{T,E}_l([\Delta \Lambda_{lm,j-1}]^{i=\infty,k} + [\delta \Lambda_{lm,j}]^{i-1,k})

.. math::
    \Delta G^{T,E}_{lm,j}=\frac{1}{g}(1+h^{T,E}_l)([\Delta \Lambda_{lm,j-1}]^{i=\infty,k} + [\delta \Lambda_{lm,j}]^{i-1,k})

Viscous components of the SLE :
-------------------------------

We define also 4 viscous component in the SLE, the ground displacement :math:`\Delta R^V_{lm}`, the geoid displacement :math:`\Delta G^V_{lm}`, the rotational ground displacement :math:`\Delta R^{T,V}_{lm}` and the rotational geoid displacement :math:`\Delta G^{T,V}_{lm}`.

.. math::
    \Delta R^V_{lm,j}=  T_l \sum_{n=0}^{j-1} h^V_l(tj-tn)\delta M_{lm,n}^{k,i}

.. math::
    \Delta G^V_{lm,j}=  T_l \sum_{n=0}^{j-1} k^V_l(tj-tn)\delta M_{lm,n}^{k,i}

.. math::
    \Delta R^{T,V}_{lm,j}=\frac{1}{g} \sum^{j-1}_{n=0} h^V_l(tj-tn)[\delta \Lambda_{lm,j}]^{i=\infty,k}

.. math::
    \Delta G^{T,V}_{lm,j}=\frac{1}{g} \sum^{j-1}_{n=0} k^V_l(tj-tn)[\delta \Lambda_{lm,j}]^{i=\infty,k}

True sediment subsidence
------------------------

.. _sed_subs:

This library was originally developed to compute effect of sediment on RSL. We considered the pure effect of sediment on RSL but also a corrected effect of sediment from water replacement. The sediment, when they are deposited, replace water and then generates an uplift induced by the diminution of ocean thickness. We choose to correct the sediment from the ocean load. 

.. math::
    \delta M_{lm,n} = \delta H_{lm,n}C_{lm,n}\rho_w

To estimate effect of sediment on RSL, you must substract the effect of the mass variation described above to the effect of sediment mass variation. 

Relative sea level variations
-----------------------------
We estimate a pure RSL where the sea level is not including variations of sediment thickness and ice thickness. 

.. math::
    \Delta SL^{i-1,k}_j = \Delta \mathcal{G}^{i-1,k}_j - \Delta R^{i-1,k}_j + \frac{\Delta \Phi ^{i-1,k}_j}{g}

The other estimation is the full RSL :

.. math::
    \Delta SL^{i-1,k}_j = \Delta \mathcal{G}^{i-1,k}_j - (\Delta R^{i-1,k}_j+\Delta H_j + \Delta I_j^{k-1}) + \frac{\Delta \Phi ^{i-1,k}_j}{g}

The resulting estimation of RSL can be compared with the ESL (only :math:`\frac{\Delta \Phi ^{i-1,k}_j}{g}`). 

Input data format
+++++++++++++++++

Mass grid format
----------------

.. _grid_format:

The different mass grid can be input as height grid, converted to mass by a simple multiplication by a defined density or as mass grid directly. The grids can be irregular or regular, they are interpolated over a sphere using stripy. These data are input as the derivative variations over time.

The topography as initial parameter is the present day topography. The initialization will update the topography according to the ice and sediment thickness. 


Implementation of Love numbers
------------------------------

.. _love:

The :cite:`dalca_2013` theory is based on the love number theory which forces us to calculate love numbers. The love numbers exits in two forms, normal mode and decay. They can also include compressible processes. We choose for computation facilities to use the love numbers computed by ALMA3 code :cite:`melini_2022`. This code is calculating incompressible decay love numbers. Benchmarking on compressible vs incompressible love numbers have demonstrated no significant difference in computed vertical displacement over 256 spherical harmonics degree. We urge you to use this code with a degree higher than 256. 

The code is working with a precise file structure for love numbers : 

| earth_model_name
| ├── h_e.dat
| ├── h_e_T.dat
| ├── h_ve.dat
| ├── h_ve_T.dat
| ├── k_e.dat
| ├── k_e_T.dat
| ├── k_ve.dat
| ├── k_ve_T.dat
| ├── l_e.dat
| ├── l_e_T.dat
| ├── l_ve.dat
| ├── l_ve_T.dat
| └── time.dat

.. table:: Corresponding love number from equation to the file names

    ======  =========================
    file    Love numbers             
    ======  =========================
    h_e     :math:`h_{\ell}^E`       
    h_e_T   :math:`h_{\ell}^{T,E}`    
    h_ve    :math:`h_{\ell}^V(t)`     
    h_ve_T  :math:`h_{\ell}^{T,V}(t)` 
    k_e     :math:`k_{\ell}^E`
    k_e_T   :math:`k_{\ell}^{T,E}`
    k_ve    :math:`k_{\ell}^V(t)`
    k_ve_T  :math:`k_{\ell}^{T,V}(t)`
    l_e     :math:`l_{\ell}^E`
    l_ve    :math:`l_{\ell}^{T,E}` 
    l_ve    :math:`l_{\ell}^V(t)` 
    l_ve_T  :math:`l_{\ell}^{T,V}(t)`
    ======  =========================

The time.dat file contains the time at which the viscous decay love numbers are computed. An example file of configurations files for ALMA3 is provided in the code supplementary files. 

.. note::
    Add the link to the ALMA3 configuration files.




