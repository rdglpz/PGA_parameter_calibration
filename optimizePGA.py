#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 10:59:23 2021

@author: Rodrigo Lopez Farias


This module contains the cost function to evaluate PGA
"""

#we load modules
import string
import random
import os
from numpy import genfromtxt


#from helper import initialize_GRASS_notebook


                         
import grass.script as gs



def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    """
    This function generates random strings.
    
    Parameters
    ----------
    size:  int
            size of the random generates string
    
    chars: array
            It contains the abecedary in uppercase
    
    
    Returns
    ________
    
    
    A random string with uppercase and digits with size
    """
    return ''.join(random.choice(chars) for _ in range(size))

def getAllPotParamas(potential):
    """
    GEt the name of the file containg the potential parameters.
    
    
    potential: string
                Name of the potential csv file
    
    Returns
    _______
    
    
    The numerical parameters potential
    """
    with open(potential) as f: s = f.readlines()[0]
    s[0:-1]
    sp = s.split(",")
    params = ""
    for w in sp[3:]:
        params = params+w+","
    return sp[2],params[:-2]





def f3(p,
       development_start, 
       development_end, 
       potential,
       devpressure,
       predictors,
       regions,
       minimize =  True,nproc = 10, repeat = 10):
    """
    Is the main cost function to evaluate the PGA parameters discount factor, compactess mean, and compactess range.
    
    
    development_start: string
                        Name of the staring developed map
    development_end: String
                        Name of the ending developed map
                        
    potential: string
        It has the name of the potential file containing the parameteres of the regression function
        
    devpressure: String
        Name of the devpressure file with parames
    
    predictors: array
        name of the map predictors or drivers that generates the suitability map
        
    regions: String
        The regions avaliable for simulating development
        
    Minimize: bookean
        Indicates if we want to maximize or minimize the cost function
        
    nproc: int
        Default processor number
        
    repeat: int
        Number of simulations per parameter combination
        
        
    Returns
    ___________
    
    
    Float value representing the optimaly of the parameters.
    
    
    """
    
    discfac,cm,cr = p
  
    
    #parmametros fijos del mapa de presion 
    gamma = 0.5
    scaling_factor = 0.1


    #resultados
    key = id_generator()
    tmpFolder = "tmp"
    if os.path.exists(tmpFolder): os.system("rm -rf "+tmpFolder)

    os.mkdir(tmpFolder)
    
    random_name = 'calib_'+key+'.csv'

    calibration_results = tmpFolder+"/"+random_name
    print(calibration_results)
    

    #ya conocemos la demanda, falta repartir la demanda de suelo entre regiones
 
    demand = "demand_2016_tr.csv"
    patches_sizes = 'patches_2011_2016'
    
    gs.run_command('r.futures.calib', 
                   development_start = development_start, 
                   development_end = development_end,
                   subregions = regions,
                   patch_sizes = patches_sizes, 
                   patch_threshold=1800, flags='s',
                   repeat = repeat, 
                   calibration_results = calibration_results, 
                   nprocs = nproc,
                   predictors=predictors,
                   devpot_params=potential, development_pressure=devpressure,
                   n_dev_neighbourhood=30, development_pressure_approach='gravity', 
                   gamma=gamma, 
                   scaling_factor=scaling_factor,
                   demand=demand, 
                   discount_factor=[discfac], 
                   compactness_mean=[cm],
                   compactness_range=[cr], 
                   num_neighbors=4, 
                   seed_search='probability', random_seed=1)
    
    info = genfromtxt(calibration_results, delimiter=',')[1:]
#    r = info[np.argmin(info[:,-1])][5]
    os.remove(calibration_results)
    minim = 1 if minimize else -1
    return info[0][-1]*minim

