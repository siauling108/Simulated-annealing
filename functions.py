# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:51:01 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import configparser

config = configparser.ConfigParser()
config.read('input.txt')

N = int( config.get('parameters', 'N') )
'''
iterations = int( config.get('parameters', 'iterations') )
nstep = int( config.get('parameters', 'nstep') )
T_min = float( config.get('parameters', 'T_min') )
T = float( config.get('parameters', 'T') )
alpha = float( config.get('parameters', 'alpha') )
'''

citynum = list(range(N))
Tem = []
accelist = []

def travel(N):
    city = np.zeros((2,N))
    for i in range(N):
        for j in (0,1):
            city[j][i] = rand.rand() #cities coordinates 
    return city

def lenght(N, city):
    global citynum
    leng=0.
    for i in range(N):
        dum=0.;
        for icoo in (0,1): #coo=coordinate, loop for each component
            iloc1=int(citynum[i])
            iloc2=int(citynum[i-1])
            dum+=(city[icoo][iloc1]-city[icoo][iloc2])**2
        leng+=np.sqrt(dum)
    return leng

def plpath(N, city):
    global citynum
    dump=np.zeros((2,N+1))
    for k in range(N):
        for i in (0,1):
            dump[i][k]=city[i][citynum[k]]
    dump[0][N]=city[0][citynum[0]]
    dump[1][N]=city[1][citynum[0]]
    plt.figure(1)
    plt.title("Followed path")
    plt.plot (dump[0],dump[1],'o-')
    plt.grid()
    plt.show()
    
def confronto():
    global accelist, Tem
    plt.title("Moves acceptance rate vs. temperature")
    plt.xlabel("Temperature (arb. units)")
    plt.ylabel("Moves acceptance rate (%)")    
    plt.plot(Tem,accelist)
    plt.legend()
    plt.grid()
    plt.show()

#--------------------------------------------------------------------------|
    
def breverse(x,y):
    global citynum, dummy
    j=max(x,y)
    i=min(x,y)
    dummy=citynum[i:j+1]
    citynum[i:j+1]=dummy[::-1] # ::-1 reverses the order of the list

def swap(i,j,N, city):
    global citynum
    i=i%N; j=j%N
    salva=citynum[i]
    citynum[i]=citynum[j]
    citynum[j]=salva

def prunegraft(x,y,z):
    global citynum, dummy
    i=min(x,y)
    j=max(x,y)
    k=z
    if k>i and k<j:
        a=citynum[0:i]
        b=citynum[j:N]
        dummy=b+a
        dummy=dummy[::-1]
        for m in range(len(dummy)):
            citynum.insert(k+m,dummy[m])
        del citynum[j+len(dummy):N+len(dummy)]
        del citynum[0:i]
    if k<=i:
        dummy=citynum[i:j]
        for m in range(len(dummy)):
            citynum.insert(k+m,dummy[m])
        del citynum[i+len(dummy):j+len(dummy)]
    if k>=j:
        dummy=citynum[i:j]
        for m in range(len(dummy)):
            citynum.insert(k+m,dummy[m])
        del citynum[i:j]
        
#--------------------------------------------------------------------------|
        
def anneal_BRev_distance(N, alpha, T, city):
    global citynum
    global Tem, accelist
    old=lenght(N, city)
    T_min = 0.1 #put them as an arg. later
    Tem=[]
    accelist=[]
    nstep=500 #number of moves 
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            breverse(x,y)
            new=lenght(N, city) 
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=lenght(N, city)
            else:
                breverse(y,x) #inverto nuovemente,
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha
        
def anneal_BRev_Metropolis(N, alpha, T, city):
    global citynum, accelist, Tem
    old=lenght(N, city)
    T_min = 0.1
    accelist=[]
    Tem=[]
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            breverse(x,y)
            new=lenght(N, city)
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=lenght(N,city)
            else:
                breverse(y,x)
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha   

        
def anneal_swap_Metropolis(N, T, alpha, city):
    global citynum, Tem, accelist
    old=lenght(N, city)
    T_min = 0.1
    Tem=[]
    accelist=[]
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            ir=int(N*rand.rand())
            swap(ir,ir+1,N, city)
            new=lenght(N, city)
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
            else:
                swap(ir+1,ir,N, city)
            
#         print "Temperatura",T
        Tem.append(T)
#         print "Prove accettate/nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha

def anneal_swap_distance(N, alpha, T, city):
    global citynum, Tem, accelist
    old=lenght(N,city)
    T_min = 0.1
    Tem=[]
    accelist=[]
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            ir=int(N*rand.rand())
            swap(ir,ir+1,N)
            new=lenght(N,city)
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
            else:
                swap(ir+1,ir,N)
        
        if acce==0: break
        Tem.append(T)
        T = T*alpha
        accelist.append(float(acce)/nstep)

def anneal_PG_Metropolis(N, alpha, T, city):
    global citynum, Tem, accelist
    old=lenght(N, city)
    T_min = 0.1
    Tem=[]
    accelist=[]
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
#si crea una citynum d'appoggio nella quale salvare la citynum originaria            
            save=[]
            for m in range(N):
                save.append(citynum[m])
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            z=int(N*rand.rand())
            prunegraft(x,y,z) #si applica il metodo
            new=lenght(N, city)
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=lenght(N, city)
            else:
                citynum=save #si torna alla citynum originaria
        
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha

def anneal_PG_distance(N, alpha, T, city):
    global citynum, Tem, accelist
    Tem=[]
    accelist=[]
    old=lenght(N, city)
    T_min = 0.1
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            save=[]
            for m in range(N):
                save.append(citynum[m])
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            z=int(N*rand.rand())
            prunegraft(x,y,z)
            new=lenght(N, city)
            
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=lenght(N, city)
            else:
                citynum=save
            
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha        