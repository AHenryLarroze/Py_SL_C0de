Usage
=====


INSTALLATION
============

.. note:: THe code is developped on windows for the moment. We want to export it on other OS but for the moment we focus on the packages. 

download 

To install the library, run in command line : 
.. code-block::

    (.venv) $ conda create -n SL_C0de 
    (.venv) $ conda activate SL_C0de
    (.venv) $ conda install pip      
    (.venv) $ pip install --index-url https://test.pypi.org/simple/ --no-deps slcode 

To install the dependencies, download the requirements.txt on https://github.com/AHenryLarroze/Py_SL_C0de and run inside the file :

.. code-block::

    (.venv) $ pip install -r requirements
    (.venv) $ conda install -c conda-forge stripy

.. note:: you need to be in the activated environment to install the packages localy

For testing the module. You can donwload the script file in https://github.com/AHenryLarroze/Py_SL_C0de. You will also need cartopy to install it. 

.. code-block::

    (.venv) $ conda install -c conda-forge cartopy