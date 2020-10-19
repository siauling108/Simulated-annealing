# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:51:01 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import configparser

config = configparser.ConfigParser()
config.read('inputs.txt')

file1 = config.get('files','distances')
file2 = config.get('files','path')
file3 = config.get('files','Tem')
file4 = config.get('files','tot_accept')

N = config.getint('parameters', 'N')
rand.seed()

def travel(N):
    '''
    Creates a random initial path made of N cities.
    Parameter:
        N: number of cities to consider.
    Returns:
        path: array containing the cities coordinates.
    '''    
    path = np.zeros((2, N))
    for i in range(N):
        for j in (0, 1):
            path[j][i] = rand.rand()
            
    return path


def length(path):
    '''
    Calculates the length of the travel.
    Parameters:
        path: coordinates of the cities (array).

    Returns:
        leng: length of the travel.
    '''
    leng = 0.
    for i in range(len(path[0])):
        for icoo in (0, 1): #coo=coordinate, loop for each component
            leng += np.sqrt((path[icoo][i]-path[icoo][i-1])**2)
        
    return leng

def Temp(T, T_min, alpha):
    '''
    Calculates and returns the temperature list.
    Parameters:
        T: maximum "temperature" of the system.
        T_min: minimum temperature.
        alpha: T scaling parameter.
    Returns:
        The temperature list.
    '''
    Tem = []
    while T >= T_min:
        Tem.append(T)
        T = T*alpha
        
    return Tem


#------------------------------------------------------------------------------|

    
def breverse(r1, r2, path):
    '''
    Executes a block reverse move.
    Parameters:
        r1, r2: two randomly generated numbers between 0 and N-1.
    Returns:
        path: updated path.
    '''
    j = max(r1,r2)
    i = min(r1,r2)
    for k in range(len(path)):
        dummy = path[k][i:j+1]
        path[k][i:j+1] = dummy[::-1] 
    
    return path
    

def swap(i, j, path):
    '''
    Executes a swap move.
    Parameters:
        i, j: two randomly generated numbers between 0 and N-1.
        N: number of cities to consider.
    Returns:
        path: updated path.
    '''
    save = path[i]
    path[i] = path[j]
    path[j] = save
    
    return path


def prunegraft(x , y, z, path):
    '''
    Executes a prune and graft move.
    Parameters:
        x, y, z: three randomly generated numbers between 0 and N-1.
    Returns:
        path: updated path.

    '''
    i = min(x, y)
    j = max(x, y)
    k = z
    if k > i and k < j:
        a = path[0:i]
        b = path[j:N]
        dummy = b+a
        dummy = dummy[::-1]
        for m in range(len(dummy)):
            path.insert(k+m, dummy[m])
        del path[j+len(dummy):N+len(dummy)]
        del path[0:i]
    if k <= i:
        dummy = path[i:j]
        for m in range(len(dummy)):
            path.insert(k+m, dummy[m])
        del path[i+len(dummy):j+len(dummy)]
    if k >= j:
        dummy = path[i:j]
        for m in range(len(dummy)):
            path.insert(k+m,dummy[m])
        del path[i:j] 
        
    return path
 
       
#--------------------------------------------------------------------------|
 
       
def anneal_BRev_distance(N, Tem, path, nstep):
    '''
    Executes the annealing procedure by using block reverse moves,
    and the minimization of the distance as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        Tem: temperatures' array.
        path: coordinates of the cities.
        nstep: number of moves for a single iteration.
    Returns:
        path: updated path.
        accelist: list containing the acceptance rate for each temperature.   
    '''
    accelist = []
    old=length(path)
    for i in range(len(Tem)):
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            r1 = int(N*rand.rand())
            r2 = int(N*rand.rand())
            while r1 > (N-1) and r2 > (N-1) : # In order to avoid IndexError
                r1 = int(N*rand.rand())
                r2 = int(N*rand.rand())
            path = breverse(r1, r2, path)
            new = length(path) 
            
            if np.exp(-(new-old)/Tem[i]) >= 1:
                old = new
                acce += 1
                new = length(path)
            else:
                path=breverse(r2, r1, path) 
                
        accelist.append(100*float(acce)/nstep)
        
    return path, accelist
 
       
def anneal_BRev_Metropolis(N, alpha, T, city, T_min, nstep, citylist):
    '''
    Executes the annealing procedure by using block reverse moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citylist: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citylist: updated citylist.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old=length(N, city, citylist)
    while T > T_min:
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            x = int(N*rand.rand())
            y = int(N*rand.rand())
            while x > (N-1) and y > (N-1) : # In order to avoid IndexError
                x = int(N*rand.rand())
                y = int(N*rand.rand())
            citylist = breverse(x, y, citylist)
            new = length(N, city, citylist)
            
            if np.exp(-(new-old)/T) > rand.rand():
                old = new
                acce += 1
                new = length(N, city, citylist)
            else:
                citylist = breverse(x, y, citylist)
                
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        T = T*alpha  
        
    return citylist, accelist, Tem

        
def anneal_swap_Metropolis(N, T, alpha, city, T_min, nstep, citylist):
    '''
    Executes the annealing procedure by using swap moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citylist: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citylist: updated citylist.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old=length(N, city, citylist)
    while T > T_min:
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            ir = int(N*rand.rand())
            ir2 = int(N*rand.rand())
            while ir > (N-1) and ir2 > (N-1) : # In order to avoid IndexError
                ir = int(N*rand.rand())
                ir2 = int(N*rand.rand())
            citylist = swap(ir, ir2, citylist)
            new = length(N, city, citylist)
            
            if np.exp(-(new-old)/T) > rand.rand():
                old = new
                acce += 1
                new = length(N, city, citylist)
            else:
                citylist = swap(ir2,ir, citylist)   
                
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        T = T*alpha
        
    return citylist, accelist, Tem


def anneal_swap_distance(N, alpha, T, city, T_min, nstep, citylist):
    '''
    Executes the annealing procedure by using swap moves,
    and the minimization of the distance as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citylist: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citylist: updated citylist.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old = length(N, city, citylist)
    while T > T_min:
        acce=0; ii=0
        while ii <= nstep:
            ii += 1
            ir = int(N*rand.rand())
            ir2 = int(N*rand.rand())
            while ir > (N-1) and ir2 > (N-1) : # In order not to get and IndexError
                ir = int(N*rand.rand())
                ir2 = int(N*rand.rand())
            citylist = swap(ir, ir2, citylist)
            new = length(N, city, citylist)
            
            if np.exp(-(new-old)/T) >= 1:
                old = new
                acce += 1
                new = length(N, city, citylist)
            else:
                citylist = swap(ir2, ir, citylist)
        
        Tem.append(T)
        accelist.append(float(acce)/nstep)
        T = T*alpha
        
    return citylist, accelist, Tem


def anneal_PG_Metropolis(N, alpha, T, city, T_min, nstep, citylist):
    '''
    Executes the annealing procedure by using prune and graft moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citylist: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citylist: updated citylist.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old=length(N, city, citylist)
    while T > T_min:
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            save = []    #needed in case of move rejection
            for m in range(N):
                save.append(citylist[m])
            x = int(N*rand.rand())
            y = int(N*rand.rand())
            z = int(N*rand.rand())
            while x > (N-1) and y > (N-1) and z > (N-1) : # In order to avoid IndexError
                x = int(N*rand.rand())
                y = int(N*rand.rand())
                z = int(N*rand.rand())
            citylist = prunegraft(x, y, z, citylist) 
            new = length(N, city, citylist)
            
            if np.exp(-(new-old)/T) > rand.rand():
                old = new
                acce += 1
                new = length(N, city, citylist)
            else:
                citylist = save
        
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        T = T*alpha
        
    return citylist, accelist, Tem


def anneal_PG_distance(N, alpha, T, city, T_min, nstep, citylist):
    '''
    Executes the annealing procedure by using prune and graft moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        N: number of the cities to consider.
        alpha: temperature scaling parameter.
        T: maximum "temperature" of the system.
        city: coordinates of the cities.
        T_min: minimum temperature.
        nstep: number of moves for a single iteration.
        citylist: list containing the number associated to each city,
                 and their visiting order (i.e. the list index number).        
    Returns:
        citylist: updated citylist.
        accelist: list containing the acceptance rate for each temperature.   
        Tem: list containing the temperatures.
    '''
    Tem = []
    accelist = []
    old=length(N, city, citylist)
    while T > T_min:
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            save = []
            for m in range(N):
                save.append(citylist[m])
            x = int(N*rand.rand())
            y = int(N*rand.rand())
            z = int(N*rand.rand())
            while x > (N-1) and y > (N-1) and z > (N-1) : # In order to avoid IndexError
                x = int(N*rand.rand())
                y = int(N*rand.rand())
                z = int(N*rand.rand())
            citylist = prunegraft(x, y, z, citylist)
            new = length(N, city, citylist)
            
            if np.exp(-(new-old)/T) >= 1:
                old = new
                acce += 1
                new = length(N, city, citylist)
            else:
                citylist = save
            
        Tem.append(T)
        accelist.append(100*float(acce)/nstep)
        T = T*alpha  
        
    return citylist, accelist, Tem
