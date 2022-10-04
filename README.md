# Patch Growing Algorithm Objective Function Calibrated With Differential Evolution and an Estimation Distribution Algorithm

Programmers:

Rodrigo López Farías
Sergio I Valdez Peña
Alberto García Robledo

This repository contains the source code of the PGA objective function implemented in the novel Estimation Distribution Algorithm and the standard Differential Evolution for the Patch Growing Algorithm Calibration.


The requirements to run this example are:

* GrassGis 7.8.6 (tested in Ubuntu)
* Installation of ```r.futures``` in GrassGis (with ```g.extension ```r.futures in GrassGis```)
* Installation of ```r.object.geometry``` (with ```g.extension r.object.geometry```)
* Database ```futures_triangle_nc```
* Modules required listed in ```packages.txt```


The Python Notebook "PGA Calibration with EDA.ipynb" and  "PGA Calibration with DE.ipynb", has examples calibrating PGA parameters (Compactness, compactess range, and scaling factor) with the Estimation distribution Optimization algorithm and Diferential Evolution.

Please be careful locating the GrassGis project and substituting the information regarding the path of GrassGis, dataset and location.


If you find useful this code please cite us:

R. Lopez-Farias, S. I. Valdez and A. Garcia-Robledo, "Parameter Calibration of the Patch Growing Algorithm for Urban Land Change Simulations," 2021 Mexican International Conference on Computer Science (ENC), 2021, pp. 1-8, doi: 10.1109/ENC53357.2021.9534789.


















