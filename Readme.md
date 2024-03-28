# SL_C0de

**SL_C0de** (/lu'make/) is a python library based on the matlab code SL_code. This module allows to calculate GIA over the world but also the sediment isostasy. It includes function to calculate radial symetric auto gravitational earth visco-elastic response. This module also includes the option to resolve the Sea level equation. This library has been benchmarked with the following Matlab code (jaustermann/SLcode/tree/master/GRL_2022_proglacial_lakes).The aim of this library is to be easy to include new sources of load and can be linked to other python library.

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
