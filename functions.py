# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:51:01 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import configparser

config = configparser.ConfigParser()
config.read('input.txt')

file1 = config.get('files','distances')
file2 = config.get('files','path')
file3 = config.get('files','Tem')
file4 = config.get('files','tot_accept')

N = config.getint('parameters', 'N')
citynum = list(range(N))

def travel(N):
    '''
    Creates a random initial path made of N cities.
    Parameter:
        N: number of cities to consider.
    Returns:
        city: array containing the cities coordinates.
    '''    
    city = np.zeros((2,N))
    for i in range(N):
        for j in (0,1):
            city[j][i] = rand.rand()
    return city


def length(N, city, citynum):
    '''
    Calculates the length of the travel.
    Parameters:
        N: number of cities to consider.
        city: coordinates of the cities.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).
    Returns:
        leng: length of the travel.
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

def Tem_length(T, T_min, alpha):
    '''
    Gives the number of the temperatures considered in the annealing process.
    Parameters:
        T: maximum "temperature" of the system.
        T_min: minimum temperature.
        alpha: T scaling parameter.
    Returns:
        The number of temperatures considered.
    '''
    count=0
    while T > T_min:
        count+=1
        T=T*alpha
    return count


def get_path(N, city, citynum):
    '''
    Calculates an array with the optimized path.
    Parameters:
        N: number of cities to consider.
        city: coordinates of the cities.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).
    '''

    path=np.zeros((2,N+1))
    for k in range(N):
        for i in (0,1):
            path[i][k]=city[i][citynum[k]]
    path[0][N]=city[0][citynum[0]]
    path[1][N]=city[1][citynum[0]]  
    
    return path


#------------------------------------------------------------------------------|

    
def breverse(x,y, citynum):
    '''
    Executes the block reverse method.
    Parameters:
        x, y: two randomly generated numbers between 0 and N-1.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).
    Returns:
        citynum: updated citynum.
    '''

    j=max(x,y)
    i=min(x,y)
    dummy=citynum[i:j+1]
    citynum[i:j+1]=dummy[::-1] 
    return citynum
    

def swap(i,j,N, citynum):
    '''
    Executes the swap method.
    Parameters:
        i, j: two randomly generated numbers between 0 and N-1.
        N: number of cities to consider.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).
    Returns:
        citynum: updated citynum.
    '''

    save=citynum[i]
    citynum[i]=citynum[j]
    citynum[j]=save
    return citynum


def prunegraft(x,y,z, citynum):
    '''
    Executes the block reverse method.
    Parameters:
        x, y, z: three randomly generated numbers between 0 and N-1.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).
    Returns:
        citynum: updated citynum.

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
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citynum: updated citynum.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old=length(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            while x>(N-1) and y>(N-1) : # In order to avoid IndexError
                x=int(N*rand.rand())
                y=int(N*rand.rand())
            citynum=breverse(x,y, citynum)
            new=length(N, city, citynum) 
            
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=length(N, city, citynum)
            else:
                citynum=breverse(y,x, citynum) 
                
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        #if acce==0: break
        T = T*alpha
    return citynum, accelist, Tem
 
       
def anneal_BRev_Metropolis(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the block reverse method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citynum: updated citynum.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old=length(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            while x>(N-1) and y>(N-1) : # In order to avoid IndexError
                x=int(N*rand.rand())
                y=int(N*rand.rand())
            citynum=breverse(x,y, citynum)
            new=length(N, city, citynum)
            
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=length(N,city, citynum)
            else:
                citynum=breverse(x,y, citynum)
                
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        #if acce==0: break
        T = T*alpha   
    return citynum, accelist, Tem

        
def anneal_swap_Metropolis(N, T, alpha, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the swap method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citynum: updated citynum.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem=[]
    accelist=[]
    old=length(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            ir=int(N*rand.rand())
            ir2=int(N*rand.rand())
            while ir>(N-1) and ir2>(N-1) : # In order to avoid IndexError
                ir=int(N*rand.rand())
                ir2=int(N*rand.rand())
            citynum = swap(ir,ir2,N, citynum)
            new=length(N, city, citynum)
            
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=length(N, city, citynum)
            else:
                citynum = swap(ir2,ir,N, citynum)   
                print('w=',np.exp(-(new-old)/T))
                
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        #if acce==0: break
        T = T*alpha
    return citynum, accelist, Tem


def anneal_swap_distance(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the swap method,
    and the minimization of the distance as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citynum: updated citynum.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem=[]
    accelist=[]
    old=length(N,city, citynum)
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
            new=length(N,city, citynum)
            
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=length(N, city, citynum)
            else:
                citynum = swap(ir2,ir,N, citynum)
        
        Tem.append(T)
        accelist.append(float(acce)/nstep)
        #if acce==0: break
        T = T*alpha
    return citynum, accelist, Tem


def anneal_PG_Metropolis(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the prune and graft method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citynum: updated citynum.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old=length(N, city, citynum)
    while T > T_min:
        rand.seed(); acce=0; ii=0
        while ii <= nstep:
            ii+=1
            save=[]    #needed in case of move rejection
            for m in range(N):
                save.append(citynum[m])
            x=int(N*rand.rand())
            y=int(N*rand.rand())
            z=int(N*rand.rand())
            while x>(N-1) and y>(N-1) and z>(N-1) : # In order to avoid IndexError
                x=int(N*rand.rand())
                y=int(N*rand.rand())
                z=int(N*rand.rand())
            citynum = prunegraft(x,y,z, citynum) 
            new=length(N, city, citynum)
            
            if np.exp(-(new-old)/T) > rand.rand():
                old=new
                acce+=1
                new=length(N, city, citynum)
            else:
                citynum=save
        
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        #if acce==0: break
        T = T*alpha
    return citynum, accelist, Tem

def anneal_PG_distance(N, alpha, T, city, T_min, nstep, citynum):
    '''
    Executes the annealing procedure by using the prune and graft method,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citynum: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citynum: updated citynum.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''

    Tem = []
    accelist = []
    old=length(N, city, citynum)
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
            while x>(N-1) and y>(N-1) and z>(N-1) : # In order to avoid IndexError
                x=int(N*rand.rand())
                y=int(N*rand.rand())
                z=int(N*rand.rand())
            citynum = prunegraft(x,y,z, citynum)
            new=length(N, city, citynum)
            
            if np.exp(-(new-old)/T) >= 1:
                old=new
                acce+=1
                new=length(N, city, citynum)
            else:
                citynum=save
            
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        #if acce==0: break
        T = T*alpha  
    return citynum, accelist, Tem
