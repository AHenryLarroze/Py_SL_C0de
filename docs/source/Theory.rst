Mathematical Theory
===================

Most of our code is based on the mathematic and computational theory from :cite:p:`dalca_2013`.

.. image:: /image/Fonctionnalite_Python.*

The SLE
-------

The Relative sea level (:math:`\Delta SL`) variations is the result of the interaction between the vertical mouvement of geoïd, :math:`\Delta Gtot` (surface of the ocean) and the vertical mouvement of the ground, :math:`\Delta Rtot` (called subsidence). We can then express the equation as :

.. math::
    \Delta SL(\Theta,\Psi,t) = \Delta Gtot(\Theta,\Psi,t) - \Delta Rtot(\Theta,\Psi,t) -\Delta H - Delta I

The geoïd variation include both, variation of the geoïd :math:`\Delta\mathcal{G}tot` surface and variation of the ocean volume. The variation of the ocean volume following a conservation of the mass, denoted :math:`\frac{\Delta \Phi}{g}`.

.. math::
    \Delta Gtot(\Theta,\Psi,t)=\Delta\mathcal{G}tot (\Theta,\Psi,t)- \frac{\Delta \Phi}{g} 

Both geoïd and ground variations (:math:`\Delta Xtot`) can be decomposed into variations (:math:`\Delta X`) due to mass redistribution and varaitions induced by earth rotation (:math:`\Delta X^T`).

.. math:: 
    \Delta Xtot (\Theta,\Psi,t)=\Delta X^{T} (\Theta,\Psi,t)+\Delta X (\Theta,\Psi,t)

The resulting SLE is : 
 
.. math:: 
    \begin{aligned}
    \Delta SL(\Theta,\Psi,t) = \Delta\mathcal{G}^{T} (\Theta,\Psi,t)+\Delta\mathcal{G} (\Theta,\Psi,t)- \frac{\Delta \Phi}{g} \\  - \Delta R^{T} (\Theta,\Psi,t)-\Delta R (\Theta,\Psi,t) - \Delta H - \Delta I
    \end{aligned}

This equation shows that variations in relative sea level are the result of the interaction of geoid and ground variations induced by mass variations, plus the effect of the earth's rotation, plus respectively the redistribution of water masses between ice, sediment and ocean and variations in the earth's surface due to sedimentary input and ice. 

Conservation of mass 
--------------------

The term :math:`\frac{\Delta \Phi}{g}` follows a conservation of mass equation based on the variation of ice (:math:`\Delta I`) and ocean volume (:math:`\Delta S`). 

.. math::

    \iint_{\Omega} \Delta I \mathrm{~d} \Omega=-\frac{\rho_{\mathrm{W}}}{\rho_{\mathrm{I}}} \iint_{\Omega} \Delta S \mathrm{~d} \Omega .

:math:`\Delta S` include three variations, the variations of the sea level, the variation of ocean volume due to ice ocean interaction and the variation of the ocean surface. These three variations are expressed as follows :

.. math:: 
    \Delta S=\Delta \mathcal{S} \mathcal{L} \cdot C+\frac{\Delta \Phi}{g} C-T_0\left[C-C_0\right]

Where :math:`T_0` is the initial ocean volume and :math:`C` is the ocean function (1 in the ocean and 0 on the continent) :

.. math::
    C= \begin{cases}1 & \text { if } Z>0 \\ 0 & \text { if } Z \leq 0\end{cases}

Injecting this expression to the conservation of mass we obtain : 

.. math::
    \begin{aligned}
    \frac{\Delta \Phi}{g}= & -\frac{1}{\mathcal{A}} \frac{\rho_{\mathrm{I}}}{\rho_{\mathrm{w}}} \iint_{\Omega} \Delta I \mathrm{~d} \Omega-\frac{1}{\mathcal{A}} \iint_{\Omega} \Delta \mathcal{S} \mathcal{L} C \mathrm{~d} \Omega \\
    & +\frac{1}{\mathcal{A}} \iint_{\Omega} T_0\left[C-C_0\right] \mathrm{d} \Omega,
    \end{aligned}

with :math:`\mathcal{A} \equiv \iint_{\Omega} C \mathrm{~d} \Omega`

Behind the ocean function the variation of topography include ice and sediment thickness. The conservation term :math:`\frac{\Delta \Phi}{g}` include then the replacement of ocean by sediment. 

Development of :math:`\Delta G` and :math:`\Delta R`
----------------------------------------------------

.. _geoid_ground_variation_theory:

To determine both :math:`\Delta G` and :math:`\Delta R`, denoted from here :math:`\Delta \chi`, :cite:`peltier_1974` and :cite:`mitrovica_1989` introduce the Green's functions that describe the response of a radial symetric self gravitating sphere. The relation include a spatial and temporal convolution between the Green functions and the Load :math:`\Delta M`.

.. math::
    \Delta \chi (\Theta,\Psi,t)=\int_{-\infty}^t \iint_{\Omega} \Delta M\left(\Theta^{\prime}, \Psi^{\prime}, t^{\prime}\right) \cdot G F\left(\gamma, t-t^{\prime}\right) \mathrm{d} \Omega^{\prime} \mathrm{d} t^{\prime}
    
Where :math:`\gamma` is :math:`cos(\gamma) = cos(\theta)cos(\theta^{\prime}) + sin(\theta)sin(\theta^{\prime})cos(\psi-\psi^{\prime})`. GF here denote the Green function.

Case of a non-rotational Earth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The GFs follows the love numbers theory :cite:`love:hal-01307751`. Our code differs from the work of :cite:`dalca_2013` by using decay love numbers where the normal mode of love numbers was used (see section xx for details about the love numbers used in this code). We use the h and k love numbers and derive for both the elastic (:math:`x_{\ell}^E`) and viscous (decay, :math:`x_{\ell}^V(t)`) part.

Here we are working on two GF, for the geoïd (:math:`\phi(\gamma,t)`) and the ground (:math:`\Gamma(\gamma,t)`) vertical motion.

.. math::
    \phi(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[\delta(t)+k_{\ell}^E \delta(t)+ k_{\ell}^V(t)\right] P_{\ell}(\cos \gamma)

.. math::
    \Gamma(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[h_{\ell}^E \delta(t)+h_{\ell}^V(t)\right] P_{\ell}(\cos \gamma)


Where :math:`a` is the Earth radius, :math:`M_e` the Earth mass, :math:`g` the gravitational constant of earth and :math:`\delta(t)` is the Dirac function. For the non-rotational part, in the convolution, GFs are applied to the Load (:math:`\Delta L (\Theta,\Psi,t)`) a pure variation of masses.

.. math::
    \Delta \chi (\Theta,\Psi,t)=\int_{-\infty}^t \iint_{\Omega} \Delta L\left(\Theta^{\prime}, \Psi^{\prime}, t^{\prime}\right) \cdot G F\left(\gamma, t-t^{\prime}\right) \mathrm{d} \Omega^{\prime} \mathrm{d} t^{\prime}

Case of a rotational Earth
^^^^^^^^^^^^^^^^^^^^^^^^^^

The effect of rotation on sea level is expressed by the perturbation of Earth's rotational vector solved by using tidal love numbers :math:`k^T` and :math:`h^T` :cite:p:`milne_1998` in the GFs, for both elastic :math:`x_{\ell}^{T,E}` and viscous :math:`x_{\ell}^{T,V}(t)`.

.. math::
    \phi^T(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[\delta(t)+k_{\ell}^{T,E} \delta(t)+ k_{\ell}^{T,V}(t)\right] P_{\ell}(\cos \gamma)

.. math::
    \Gamma^T(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[h_{\ell}^{T,E} \delta(t)+h_{\ell}^{T,V}(t)\right] P_{\ell}(\cos \gamma)

Where 

For the rotational Earth convolution a rotational potential is defined as :math:`\Lambda(\Theta,\Psi,t_j)`. The equations behind are described in :cite:`milne_1998` and are not developed here.

Resolution of temporal and spatial convolution
----------------------------------------------

Spatial convolution
^^^^^^^^^^^^^^^^^^^

The spatial convolution is resolved using the spherical harmonic transformation. For a function :math:`\chi(\Theta,\Psi,t)`, we can define spherical harmonics coefficients :math:`\chi_{lm}(t)`, where :math:`l` is the degree and :math:`m` is the order of the associated Legendre polynomial (:math:`Y_{lm}(\Theta,\Psi)`) :

.. math::
    \mathcal{X} (\Theta,\Psi,t)=\sum_{lm} \mathcal{X} _{lm}(t)Y_{lm}(\Theta,\Psi)

with :math:`\sum_{lm}=\sum_{l=0}^{\infty} \sum_{m=-l}^{m=l}`, for a degree l there is 2l+1 order.  

In the spectral domain the convolution can be solved : 

.. math::
    \iint _{\Omega} \sum_{l=0}^{infty} \mathcal{X} (\Theta',\Psi',t)P_l(cos\gamma') \,d \Omega = T_l \sum_{lm} \mathcal{X} _{lm} (t) Y_{lm}(\Theta,\Psi)
    
.. _T_definition:

With :math:`T_l = \frac{4\pi a^2}{2l+1}`


Temporal convolution
^^^^^^^^^^^^^^^^^^^^

The resolution of temporal convolution is performed by a Heaviside distribution of the load :math:`\mathcal{H} (t)`. 

 .. math:: 
    \mathcal{H} (t) = \left\{
        \begin{array}{ll}
            1 & \mbox{si t>0} \\
            \varnothing   & \mbox{si t=0} \\
            0 & \mbox{si t<0} \\
        \end{array}
    \right.

The Heavyside distributed load is : 

.. math::
    \Delta L(\Theta,\Psi,t)=\sum_{n=0}^{N} \delta L(\Theta,\Psi,t_n)\mathcal{H} (t-t_n)

Resolution of the convolutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Applying the temporal convolution resolution :

.. math::
    \Delta \chi= \iint_{\Omega} \sum_{n=0}^{N} \delta M(\Theta,\Psi,t_n) \int_{-\infty}^t  \mathcal{H} (t-t_n) \cdot G F\left(\gamma, t-t^{\prime}\right) \mathrm{d} \Omega^{\prime} \mathrm{d} t^{\prime}

and : 

.. math::
    \int_{-\infty}^t  \mathcal{H} (t-t_n) \cdot G F\left(\gamma, t-t^{\prime}\right) \mathrm{d} t^{\prime} = I G F\left(\gamma, t-t_n \right)

with IGF the time integration of the GF.

We have then : 

.. math::
    \Delta \chi (\Theta,\Psi,t)= \iint_{\Omega} \sum_{n=0}^{j-1} \delta M(\Theta,\Psi,t_n) \cdot I G F\left(\gamma, t_j -t_n \right)

By application of the spatial convolution solution : 

.. math::
     \Delta \chi (\Theta,\Psi,t)= \sum_{lm} T_l \sum_{n=0}^{j-1} \delta M_{lm} (t_n) Y_{lm}(\Theta,\Psi) \cdot I G F\left(\gamma, t-t_n \right)

The respective IGF are :

.. math::
    I \phi(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[1+k_{\ell}^E+ K_{\ell}^V(t)\right]

.. math::
    I \Gamma(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[h_{\ell}^E +H_{\ell}^V(t)\right]

.. math::
    I \phi^T(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[1+k_{\ell}^{T,E} + K_{\ell}^{T,V}(t)\right]

.. math::
    I \Gamma^T(\gamma, t)=\frac{a g}{M_{\mathrm{e}}} \sum_{\ell=0}^{\infty}\left[h_{\ell}^{T,E}+H_{\ell}^{T,V}(t)\right]

Where K and H are the integrated love numbers between 0 and the considered time t. 

The resulting SLE is :

.. math::
    \Delta \mathcal{SL}(\Theta,\Psi,t) = \int_{-\infty}^{t_j} \iint_\Omega \Delta L (\Theta ^\prime, \Psi ^\prime,t ^\prime) \cdot [\frac{\Phi(\gamma, t - t^\prime)}{g} - \Gamma(\gamma,t-t^\prime)] \mathrm{d} \Omega ^\prime \mathrm{d} t^\prime - \Delta H(\Theta,\Psi,t) - \Delta I(\Theta,\Psi,t)


