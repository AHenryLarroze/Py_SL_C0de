# SL_C0de

SL_C0de is a python library based on the matlab code SL_code. This module is used to calculate GIA over the world but alse the sediment isostasy. It include function to calculate radial symetric auto gravitational earth visco-elastic response. This module also include the possibility to resolve the Sea level equation based on the publication of Mitrovica and Austermann. The aim of this library is to be easy to include new sources of load and link to other python library.

## install : 

### install the main package : 

```
    $ conda create -n SL_C0de 
    $ conda activate SL_C0de
    $ conda install pip      
    $ pip install -i https://test.pypi.org/simple/ --no-deps slcode==0.4.1
```

### install requirements : 

```
    $ pip install -r requirements
    $ conda install cartopy
```

## Easy use : 

Download the /Data file and the /script/SL_C0de_run.ipynb

require panda