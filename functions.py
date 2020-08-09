# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:51:01 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt


def lung():
    global N,city,lista
    leng=0.
    for i in range(N):
        dum=0.;
        for icoo in (0,1): #coo=coordinate, loop for each component
            iloc1=int(lista[i])
            iloc2=int(lista[i-1])
            dum+=(city[icoo][iloc1]-city[icoo][iloc2])**2
        leng+=np.sqrt(dum)
    return leng

def plpath():
    global city,N,lista
    dump=np.zeros((2,N+1))
    for k in range(N):
        for i in (0,1):
            dump[i][k]=city[i][lista[k]]
    dump[0][N]=city[0][lista[0]]
    dump[1][N]=city[1][lista[0]]
    plt.figure(1)
    plt.title("Percorso tra le citta'")
    plt.plot (dump[0],dump[1],'o-')
    plt.grid()
    plt.show()
    
def confronto():
    global accelist, Tem
    plt.title("Mosse accettate per numero di passi in funzione della temperatura")
    plt.xlabel("Temperatura(unita' arbitrarie)")
    plt.ylabel("Mosse accettate (%)")    
    plt.plot(Tem,accelist)
    plt.legend()
    plt.grid()
    plt.show()

#--------------------------------------------------------------------------|
    
def breverse(x,y):
    global lista, dummy
    j=max(x,y)
    i=min(x,y)
    dummy=lista[i:j+1]
    lista[i:j+1]=dummy[::-1] # ::-1 reverses the order of the list

def swap(i,j):
    global N, city, lista
    i=i%N; j=j%N
    salva=lista[i]
    lista[i]=lista[j]
    lista[j]=salva

def prunegraft(x,y,z):
    global lista, dummy
    i=min(x,y)
    j=max(x,y)
    k=z
    if k>i and k<j:
        a=lista[0:i]
        b=lista[j:N]
        dummy=b+a
        dummy=dummy[::-1]
        for m in range(len(dummy)):
            lista.insert(k+m,dummy[m])
        del lista[j+len(dummy):N+len(dummy)]
        del lista[0:i]
    if k<=i:
        dummy=lista[i:j]
        for m in range(len(dummy)):
            lista.insert(k+m,dummy[m])
        del lista[i+len(dummy):j+len(dummy)]
    if k>=j:
        dummy=lista[i:j]
        for m in range(len(dummy)):
            lista.insert(k+m,dummy[m])
        del lista[i:j]
        
#--------------------------------------------------------------------------|
        
def anneal_BRev_distance():
    global N, city, lista, alpha, T, Tem, accelist
    old=lung()
    T_min = 0.1
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
            new=lung() 
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=lung()
            else:
                breverse(y,x) #inverto nuovemente,
            pat.append(new)     #tornando al caso iniziale
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha
        
def anneal_BRev_Metropolis():
    global N, city, lista, alpha, T, accelist, Tem
    old=lung()
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
            new=lung()
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=lung()
            else:
                breverse(y,x)
            pat.append(new)
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha   

        
def anneal_swap_Metropolis():
    global N, city, lista, alpha, T, Tem, accelist, acce
    old=lung()
    T_min = 0.1
    Tem=[]
    accelist=[]
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            ir=int(N*rand.rand())
            swap(ir,ir+1)
            new=lung()
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
            else:
                swap(ir+1,ir)
            pat.append(new)
#         print "Temperatura",T
        Tem.append(T)
#         print "Prove accettate/nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha

def anneal_swap_distance():
    global N, city, lista, alpha, T, Tem, accelist
    old=lung()
    T_min = 0.1
    Tem=[]
    accelist=[]
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            ir=int(N*rand.rand())
            swap(ir,ir+1)
            new=lung()
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
            else:
                swap(ir+1,ir)
            pat.append(new)
        if acce==0: break
        Tem.append(T)
        T = T*alpha
        accelist.append(float(acce)/nstep)

def anneal_PG_Metropolis():
    global N, city, lista, alpha, T, Tem, accelist
    old=lung()
    T_min = 0.1
    Tem=[]
    accelist=[]
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
#si crea una lista d'appoggio nella quale salvare la lista originaria            
            save=[]
            for m in range(N):
                save.append(lista[m])
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            z=int(N*rand.rand())
            prunegraft(x,y,z) #si applica il metodo
            new=lung()
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=lung()
            else:
                lista=save #si torna alla lista originaria
            pat.append(new)
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha

def anneal_PG_distance():
    global N, city, lista, alpha, T, save, Tem, accelist
    Tem=[]
    accelist=[]
    old=lung()
    T_min = 0.1#0.00001
    nstep=500
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            save=[]
            for m in range(N):
                save.append(lista[m])
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            z=int(N*rand.rand())
            prunegraft(x,y,z)
            new=lung()
            
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=lung()
            else:
                lista=save
            
            pat.append(new)
#         print "T",T
        Tem.append(T)
#         print "Acce/Nstep", float(acce)/nstep
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha        