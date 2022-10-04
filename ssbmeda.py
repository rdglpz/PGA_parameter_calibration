##############################################################################
####S. Ivvan Valdez
####November 20 2019
####Implementation of the Estimation of Distribution Algorithm for expensive objective functions
####This file is a library with the main function for the EDA-EOF
####The algorithm is based on the on the Normal selection
####It considers a maximization case

from numpy.random import seed
from numpy.random import rand
from numpy.random import randn
import numpy.random as npr
from scipy.stats import spearmanr as spr
import numpy as np
from time import time

import importlib
#import tpc as problems
import empsel as esel

print("----2")



# importlib.reload(problems)
# importlib.reload(esel)
# reload(problems)
# reload(esel)
seed(int(time()))

def popIni(dim,size):
    X=rand(size,dim)
    return X

def evaluate(X, problem,minmax=1.0):
    nrow=len(X)
    F=np.zeros(nrow)
    i=0
    for  x in X:
        F[i]=minmax*problem(x)
        i+=1
    return F

def genNInd3(X,F,xelite,p,ibest,N,narch,problem,minmax,u):

    if 'counter2' not in  genNInd3.__dict__ or set==1:
        genNInd3.counter2=0
    if 'normsd1' not in  genNInd3.__dict__ or set==1:
        genNInd3.normsd1=-1
    if 'normsd2' not in  genNInd3.__dict__ or set==1:
        genNInd3.normsd2=-1
    if 'felite' not in  genNInd3.__dict__ or set==1:
        genNInd3.felite=-1e100

#INTENTARE MANTENER LA ULTIMA NORMA EXITOSA DE LA DESVIACION ESTANDARD
    success=False
    if genNInd3.felite<F[ibest[0]] and genNInd3.normsd1!=-1:
        success=True
    genNInd3.felite=F[ibest[0]]
    count=10
    nvar=np.shape(X)
    npop=nvar[0]
    nvar=nvar[1]
    ntrials=2000 #Funciono mejor con 1500
    nb=-1
    xn=np.zeros((N,nvar))
    meanX=sum(np.transpose(np.transpose(X)*p))
    sdX=np.sqrt(sum(np.transpose(np.transpose((X-meanX)**2)*p)))

    genNInd3.counter2+=1
    if success==True:
       genNInd3.counter2-=1

    increment=True
    # if  genNInd3.counter2==count or genNInd3.counter2==0:
    if  genNInd3.counter2>=(int(count/2)+1) or genNInd3.counter2==0:
        mult=1.05 #Ahorita lo voy a bajar
        sdX=sdX*mult
        scale1=np.max((genNInd3.normsd1/np.linalg.norm(sdX),1.0))
        print('Increment',success,'++++++')
        if  genNInd3.counter2==count:
            genNInd3.counter2=0
    else:
        mult=0.6
        sdX=sdX*mult
        scale1=np.min((genNInd3.normsd1/np.linalg.norm(sdX),1.0))
        increment=False
        print('Decrement',success,'------')


##Repetir el caso anterior si le servia!!! con 0.75 y 1.05 de expansion y reduccion mejores resultados hasta Ahorita
##Llega a lo mas bajo en un 40% mas o menos pero mas bajo que los casos anteriores [-0.01473454, -0.02773657, -3.56831343, -0.05539519, -0.03063693
##Para 1.05 y 0.9 [-0.05022853, -3.54329088, -0.03647488, -3.54376348, -0.02462781,-6.17042122, -6.17040435, -0.04200643, -0.05062388, -3.54263786]
#Para 1.05 y 0.7 [-0.05612015, -3.54409543, -0.03888253,
    genNInd3.normsd1=np.linalg.norm(sdX*scale1)


    for i in range(N):
        xn[i,]=X[ibest[0],]+randn(nvar)*sdX

    # Fantes=np.zeros(N)

    uf=u #np.zeros(nvar)
    ut=np.dot(X,u)
    crmax=spr(ut,F)[0]
    for i in range(ntrials):
        ut=uf+0.25*randn(nvar)/np.sqrt(nvar)
        ut=ut/np.linalg.norm(ut)
        u0=np.dot(X,ut)
        cr=spr(u0,F)[0]
        if abs(cr)>abs(crmax):
            crmax=cr
            uf=ut
            if cr<0:
                uf=-u
                crmax=-cr
        # if abs(crmax)>0.95:
        #     break
    #Mejor resultado para 0.51 hasta ahora [-0.01473454, -0.02773657, -3.56831343, -0.05539519, -0.03063693
    #Con 0.61 sale mucho peor [-0.03168507, -6.1714098 , -3.54252487, -3.54271004, -3.54532564,-3.5432036 , -6.17041551,
    #Usando el anterior u anterior[-0.03302358, -0.05278755, -0.05488436, -6.17789788, -3.54385873
    #Usando una combinacion de los vectores anteriores -6.17038224e+00, -8.28068426e-03, -5.03804023e-03, -2.05673455e-03,-6.17038203e+00, -5.63562625e-03, -3.54094873e+00, -8.42778139e-03,-6.31097702e-03, -6.17038208e+00]
    if abs(crmax)>-1: #0.51:
        for i in range(N):
            xn[i,]=X[ibest[0],]+randn(nvar)*sdX*scale1
        # for i in range(N):
        #     Fantes[i]=minmax*problem(xn[i,])
        uf=uf/np.linalg.norm(uf)
        u=u/np.linalg.norm(u)
        uf=0.05*uf+0.95*u
        uf=uf/np.linalg.norm(uf)
        st=np.dot(X[ibest[0]]-X[ibest[narch-1]],uf)
        fs=0.6-2.0*(crmax-1.0)**2 #Mejor hasta ahora
        scale2=1.0;
        if increment==True:
            scale2=np.max((genNInd3.normsd2/st,1.0))
        else:
            scale2=np.min((genNInd3.normsd2/st,1.0))
        genNInd3.normsd2=st*scale2
        for i in range(N):
            xn[i,]+=(fs*randn()+fs)*st*uf
            xn[i,xn[i,]<0.0]=0.0
            xn[i,xn[i,]>1.0]=1.0


    return np.reshape(xn,(N,nvar)),nb,uf


#####------------------------------------#####
####Default algorithm parameters
def EDAEOF(tproblem,dim,npop=-1,ngen=-1,narch=-1,selection=esel.ebintour,pressure=1.0,printing=1,maxgen=100000,minVar=1.0e-5,maxEval=5e4,maxFval=0,runfile="test0001.out",outfile="exec0001.dat",appendfile="aexec0001.dat"):
    if npop==-1:
        npop=2*dim
    if ngen==-1:
        ngen=dim
    if narch==-1:
        narch=5*npop


    cvar=1e100
    ite=0
    minmax=-1.0
    elitec=0
    resets=0;
    ####Algorithm
    X=popIni(dim,npop)
    F=evaluate(X,tproblem,minmax=minmax)
    ibest,ep=selection(F,set=1,pressure=pressure)
    xelite=X[ibest[0]]
    felite=F[ibest[0]]
    freset=felite
    neval=len(F)
    ite+=1
    # nbacu=0
    # nbenter=0
    u=npr.rand(dim)
    outfile=open(outfile,'a')
    runfile=open(runfile,'a')
    appendfile=open(appendfile,'a')


#    print >>runfile,"Gen=",ite,"Fbest=",F[ibest[0]],"cstd=",cvar,"neval=",neval,"elitec=",elitec,"resets=",resets
    print("Gen=",ite,"Fbest=",F[ibest[0]],"cstd=",cvar,"neval=",neval,"elitec=",elitec,"resets=",resets,file=runfile)
    #print >>runfile,'X=',X
    print('X=',X,file=runfile)
    #print >>runfile,'F=',F
    print('F=',F,file=runfile)
    #print >>outfile, ite*np.ones(len(F)),neval*np.ones(len(F)),F[ibest],X[ibest]
    print(ite*np.ones(len(F)),neval*np.ones(len(F)),F[ibest],X[ibest],file=outfile) 
    
    while  ite<maxgen and cvar>minVar and neval<maxEval and F[ibest[0]]<maxFval:
        print("iteration = ",ite)
        xt,nb,u=genNInd3(X,F,xelite,ep,ibest,ngen,narch,tproblem,minmax=minmax,u=u)
        ft=evaluate(xt,tproblem,minmax=minmax)
        print("aptval:", ft)
        neval+=len(ft)
        if (abs(felite-F[ibest[0]])>abs(F[0])*1e-4):
            elitec=0
            xelite=X[ibest[0]]
            felite=F[ibest[0]]
        else:
            elitec+=1

        if elitec==int(resets*10+15):
            resets+=1
            xelite=X[ibest[0:resets]]
            felite=F[ibest[0:resets]]
            X=popIni(dim,npop-1)
            X=X-0.5
            X=(1-0.15*resets)*X
            for i in range(dim):
                X[:,i]+=xelite[0,i]
                X[X[:,i]<0.0,i]=0.0
                X[X[:,i]>1.0,i]=1.0
            F=evaluate(X,tproblem,minmax=minmax)
            X=np.concatenate((X, np.reshape(xelite,(resets,dim))))
            # F=np.concatenate(F,np.array([felite])
            F=np.concatenate((F,felite),axis=None)
            ibest,ep=selection(F,set=1,pressure=pressure)
            u=npr.rand(dim)
            xelite=X[ibest[0]]
            felite=F[ibest[0]]
            neval+=(len(F)-1)
            elitec=0
            if (resets==narch or felite==freset):
                break
            else:
                freset=felite
            #print >>runfile, 'Reset!!!',"elitec=",elitec,"resets=",resets,"var=",np.var(X)
            print('Reset!!!',"elitec=",elitec,"resets=",resets,"var=",np.var(X),file=runfile)
        X=X[ibest]
        F=F[ibest]
        if len(ibest)>narch:
            F=F[0:narch]
            X=X[0:narch]
        X=np.concatenate((X,xt))
        F=np.concatenate((F,ft))

        ibest,ep=selection(F,pressure=pressure)
        #print >>outfile, ite*np.ones(len(F)),neval*np.ones(len(F)),F[ibest],X[ibest]
        print(ite*np.ones(len(F)),neval*np.ones(len(F)),F[ibest],X[ibest],file=outfile)
        cvar=np.std(F[ep>0])
        if printing>=2:
            #print >>runfile,
            print("Gen=",ite,"Fbest=",F[ibest[0]],"cstd=",cvar,"neval=",neval,"elitec=",elitec,"resets=",resets,file=runfile)
            #print >>runfile,'X=',X
            print('X=',X,file=runfile)
            #print >>runfile,'F=',F
            print('F=',F,file=runfile)
        ite+=1
    if printing>=1:
        #print >>runfile,"Gen=",ite,"Fbest=",F[ibest[0]],"cstd=",cvar,"neval=",neval,"Xbest=",X[ibest[0]]
        print("Gen=",ite,"Fbest=",F[ibest[0]],"cstd=",cvar,"neval=",neval,"Xbest=",X[ibest[0]],file=runfile)
        #print >>appendfile,"Ite=",ite,"Fbest=",F[ibest[0]],"Xbest=",X[ibest[0]],"cstd=",cvar,"neval=",neval
        print("Ite=",ite,"Fbest=",F[ibest[0]],"Xbest=",X[ibest[0]],"cstd=",cvar,"neval=",neval,file=appendfile)
        #print >>outfile, ite*np.ones(len(F)),neval*np.ones(len(F)),F[ibest],X[ibest]
        print(ite*np.ones(len(F)),neval*np.ones(len(F)),F[ibest],X[ibest],file=outfile)

    outfile.close()
    runfile.close()
    appendfile.close()

    return X[ibest],F[ibest],neval,resets
