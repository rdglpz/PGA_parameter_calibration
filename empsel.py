##############################################################################
####S. Ivvan Valdez
####November 20 2019
####Implementation of some empirical selection distribution (ESD) functions
####The functions computes the probability of an individual to be selected.
####The esample function selects n individuals given the ESD

import numpy as np
from scipy.stats import norm
import copy

#import matplotlib.pyplot as plt
##-------IMPORTANT COMMENTS!!!!
#THE SELECTION METHODS ARE MADE FOR MAXIMIZATION

def ebintour(F,set=0,pressure=1.0):
    p=np.zeros(len(F))
    ibest=np.argsort(F)[::-1]
    p[ibest]=np.arange(start=len(F),stop=0,step=-1)
    p=p**pressure
    p=p/sum(p)
    return ibest,p

def ebintrunc(F,set=0,pressure=1.0):
    if 'threshold' not in  ebintrunc.__dict__ or set==1:
        ebintrunc.threshold=min(F)
    p=np.zeros(len(F))
    ibest=np.argsort(F)[::-1]
    p[ibest]=np.arange(start=len(F),stop=0,step=-1)
    mp=min(p[F>=ebintrunc.threshold])
    p=p-mp+1.0
    if sum(p>0)>4:
        p[p<=0]=0
    else:
        p=p+mp-1.0
    p=p**pressure
    p=p/sum(p)
    ebintrunc.threshold=min(F[p>0])
    # print(ebintrunc.threshold,p)
    return ibest,p

def etrunc(F,set=0,pressure=1.0):
    if 'threshold' not in  etrunc.__dict__ or set==1:
        etrunc.threshold=min(F)
    p=np.zeros(len(F))
    ibest=np.argsort(F)[::-1]
    npop=len(F)
    nsel=int(npop/2)
    p[ibest[0:nsel]]=1.0/nsel
    p=p/sum(p)
    etrunc.threshold=min(F[p>0])
    return ibest,p


def ebintrunc2(F,set=0,pressure=1.0):
    if 'threshold' not in  ebintrunc2.__dict__ or set==1:
        ebintrunc2.threshold=min(F)
    p=np.zeros(len(F))
    ibest=np.argsort(F)[::-1]
    p[ibest]=np.arange(start=len(F),stop=0,step=-1)
    t2=F[ibest[int(len(F)/2)]]
    ebintrunc2.threshold=max(ebintrunc2.threshold,t2)
    mp=min(p[(F>=ebintrunc2.threshold)])
    p=p-mp+1.0
    if sum(p>0)>4:
        p[p<=0]=0
    else:
        p=p+mp-1.0
    p=p**pressure
    p=p/sum(p)
    ebintrunc2.threshold=min(F[p>0])
    print("threshold=",ebintrunc2.threshold,"t=",t2)
    return ibest,p

def ebintrunc3(F,set=0,pressure=1):
    p=np.zeros(len(F))
    ibest=np.argsort(F)[::-1]
    p[ibest]=np.arange(start=len(F),stop=0,step=-1)
    ebintrunc3.threshold=F[ibest[int(pressure-1)]]
    mp=min(p[(F>=ebintrunc3.threshold)])
    p=p-mp+1.0
    if sum(p>0)>4:
        p[p<=0]=0
    else:
        p=p+mp-1.0
    p=p**pressure
    p=p/sum(p)
    print("threshold=",ebintrunc3.threshold)
    return ibest,p


def enorm(F, set=0,pressure=2.5):
    p=(norm(0, 1).cdf(np.arange(0,len(F)+1,1)*pressure/len(F))-0.5)
    pt=p/p[-1]
    p=pt[1:(len(F)+1)]-pt[0:len(F)]
    p=p/sum(p)
    pt=copy.deepcopy(p)
    ibest=np.argsort(F)[::-1]
    p[ibest]=pt
    return ibest,p

def esample(size,p,replace=True):
    s=np.random.choice(np.arange(0,len(p),1),size=size,replace=replace,p=p)
    return p
