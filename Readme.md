# :math:`SL_{C0de}`'

:math:`SL_{C0de}`' is a python library based on the matlab code SL_code. This module is used to calculate GIA over the world but alse the sediment isostasy. It include function to calculate radial symetric auto gravitational earth visco-elastic response. This module also include the possibility to resolve the Sea level equation based on the publication of Mitrovica and Austermann. The aim of this library is to be easy to include new sources of load and link to other python library.

## install : 

### install the main package : 

.. code-block::
    
    (.venv) $ conda create -n SL_C0de 
    (.venv) $ conda activate SL_C0de
    (.venv) $ conda install pip      
    (.venv) $ pip install -i https://test.pypi.org/simple/ --no-deps slcode==0.4.1

### install requirements : 

.. code-block::

    (.venv) $ pip install -r requirements
    (.venv) $ conda install cartopy
	