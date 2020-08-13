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

N = config.getint('parameters', 'N')
citynum = list(range(N))

def travel(N):
    '''
    Creates a random initial path made of N cities.
    Parameter:
        N: number of cities to consider.
    '''    
    city = np.zeros((2,N))
    for i in range(N):
        for j in (0,1):
            city[j][i] = rand.rand() #cities coordinates 
    return city


def lenght(N, city, citynum):
    '''
    Calculates the lenght of the travel.
    Parameters:
        N: number of cities to consider.
        city: coordinates of the cities.
    '''

    leng=0.
    for i in range(N):
        dum=0.;
        for icoo in (0,1): #coo=coordinate, loop for each component
            iloc1=int(citynum[i])
            iloc2=int(citynum[i-1])
            dum+=(city[icoo][iloc1]-city[icoo][iloc2])**2
        leng+=np.sqrt(dum)
    return leng


def plpath(N, city, citynum):
    '''
    Plots the followed path.
    Parameters:
        N: number of cities to consider.
        city: coordinates of the cities.
    '''

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
 
    
def acceptance_plot(accelist, Tem):
    '''
    Plots the acceptance rate as a function of the temperature.
    '''

    plt.title("Moves acceptance rate vs. temperature, 1st iteration")
    plt.xlabel("Temperature (arb. units)")
    plt.ylabel("Moves acceptance rate (%)")    
    plt.plot(Tem,accelist)
    plt.legend()
    plt.grid()
    plt.show()

#--------------------------------------------------------------------------|
    
def breverse(x,y, citynum):
    '''
    Executes the block reverse method.
    Parameters:
        x, y: two randomly generated numbers between 0 and N.
    '''

    j=max(x,y)
    i=min(x,y)
    dummy=citynum[i:j+1]
    citynum[i:j+1]=dummy[::-1] # ::-1 reverses the order of the list
    return citynum
    

def swap(i,j,N, citynum):
    '''
    Executes the swap method.
    Parameters:
        i: randomly generated number between 0 and N.
        j: i+1.
    '''

    #i=i%N; j=j%N
    save=citynum[i]
    citynum[i]=citynum[j]
    citynum[j]=save
    return citynum


def prunegraft(x,y,z, citynum):
    '''
    Executes the block reverse method.
    Parameters:
        x, y, z: three randomly generated numbers between 0 and N.
    '''

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
    return citynum
 
       
#--------------------------------------------------------------------------|
        
def anneal_BRev_distance(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the block reverse method,
    and the minimization of the distance as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: rescaling parameter for the temperature.
        T: "temperature" of the system.
        city: coordinates of the cities.
        T_min: temperature minimum.
        nstep: number of moves for a single iteration.
        Tem: list containing the temperatures.
        accelist: list containing the acceptance rate for each temperature.        
    '''
    Tem = []
    accelist = []
    old=lenght(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            citynum=breverse(x,y, citynum)
            new=lenght(N, city, citynum) 
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=lenght(N, city, citynum)
            else:
                citynum=breverse(y,x, citynum) #inverto nuovamente,
                
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha
    return citynum, accelist, Tem
 
       
def anneal_BRev_Metropolis(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the block reverse method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: rescaling parameter for the temperature.
        T: "temperature" of the system.
        city: coordinates of the cities.
        T_min: temperature minimum.
        nstep: number of moves for a single iteration.
        Tem: list containing the temperatures.
        accelist: list containing the acceptance rate for each temperature.        
    '''
    Tem = []
    accelist = []
    old=lenght(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            citynum=breverse(x,y, citynum)
            new=lenght(N, city, citynum)
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=lenght(N,city, citynum)
            else:
                citynum=breverse(x,y, citynum)
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha   
    return citynum, accelist, Tem

        
def anneal_swap_Metropolis(N, T, alpha, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the swap method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: rescaling parameter for the temperature.
        T: "temperature" of the system.
        city: coordinates of the cities.
        T_min: temperature minimum.
        nstep: number of moves for a single iteration.
        Tem: list containing the temperatures.
        accelist: list containing the acceptance rate for each temperature.        
    '''
    Tem=[]
    accelist=[]
    old=lenght(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            ir=int(N*rand.rand())
            ir2=int(N*rand.rand())
            while ir>(N-1) and ir2>(N-1) : # In order not to get and IndexError
                ir=int(N*rand.rand())
                ir2=int(N*rand.rand())
            citynum = swap(ir,ir2,N, citynum)
            new=lenght(N, city, citynum)
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=lenght(N, city, citynum)
            else:
                citynum = swap(ir2,ir,N, citynum)            
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha
    return citynum, accelist, Tem


def anneal_swap_distance(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the swap method,
    and the minimization of the distance as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: rescaling parameter for the temperature.
        T: "temperature" of the system.
        city: coordinates of the cities.
        T_min: temperature minimum.
        nstep: number of moves for a single iteration.
        Tem: list containing the temperatures.
        accelist: list containing the acceptance rate for each temperature.        
    '''
    Tem=[]
    accelist=[]
    old=lenght(N,city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            ir=int(N*rand.rand())
            ir2=int(N*rand.rand())
            while ir>(N-1) and ir2>(N-1) : # In order not to get and IndexError
                ir=int(N*rand.rand())
                ir2=int(N*rand.rand())
            citynum = swap(ir,ir2,N, citynum)
            new=lenght(N,city, citynum)
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=lenght(N, city, citynum)
            else:
                citynum = swap(ir2,ir,N, citynum)
        
        if acce==0: break
        Tem.append(T)
        T = T*alpha
        accelist.append(float(acce)/nstep)
    return citynum, accelist, Tem


def anneal_PG_Metropolis(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the prune and graft method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: rescaling parameter for the temperature.
        T: "temperature" of the system.
        city: coordinates of the cities.
        T_min: temperature minimum.
        nstep: number of moves for a single iteration.
        Tem: list containing the temperatures.
        accelist: list containing the acceptance rate for each temperature. 
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).
    '''
    Tem = []
    accelist = []
    old=lenght(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
#si crea una lista d'appoggio nella quale salvare la citynum originaria            
            save=[]
            for m in range(N):
                save.append(citynum[m])
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            z=int(N*rand.rand())
            citynum = prunegraft(x,y,z, citynum) #si applica il metodo
            new=lenght(N, city, citynum)
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=lenght(N, city, citynum)
            else:
                citynum=save #si torna alla citynum originaria
        
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha
    return citynum, accelist, Tem

def anneal_PG_distance(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the prune and graft method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: rescaling parameter for the temperature.
        T: "temperature" of the system.
        city: coordinates of the cities.
        T_min: temperature minimum.
        nstep: number of moves for a single iteration.
        Tem: list containing the temperatures.
        accelist: list containing the acceptance rate for each temperature. 
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number)
    '''

    Tem = []
    accelist = []
    old=lenght(N, city, citynum)
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
            citynum = prunegraft(x,y,z, citynum)
            new=lenght(N, city, citynum)
            
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=lenght(N, city, citynum)
            else:
                citynum=save
            
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        if acce==0: break
        T = T*alpha  
    return citynum, accelist, Tem
