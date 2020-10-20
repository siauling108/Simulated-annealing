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
        path: coordinates of the cities (array).
    Returns:
        path: updated path.
    '''
    j = max(r1,r2)
    i = min(r1,r2)
    for k in range(len(path)):
        dummy = path[k][i:j+1]
        path[k][i:j+1] = dummy[::-1] 
    
    return path
    

def swap(r1, r2, path):
    '''
    Executes a swap move.
    Parameters:
        r1, r2: two randomly generated numbers between 0 and N-1.
        path: coordinates of the cities (array).
    Returns:
        path: updated path.
    '''
    for i in range(len(path)):
        save = path[i][r1]
        path[i][r1] = path[i][r2]
        path[i][r2] = save
    
    return path


def prunegraft(r1 , r2, r3, path):
    '''
    Executes a prune and graft move.
    Parameters:
        r1 , r2, r3: three randomly generated numbers between 0 and N-1.
        path: coordinates of the cities (array).
    Returns:
        path: updated path.

    '''
    i = min(r1, r2)
    j = max(r1, r2)
    k = r3
    
    for d in range(len(path)):
        path_l = path[d].tolist()
        if k > i and k < j:
            a = path_l[0:i]
            b = path_l[j:N]
            dummy = b+a
            dummy = dummy[::-1]
            for m in range(len(dummy)):
                path_l.insert(k+m, dummy[m])
            del path_l[j+len(dummy):N+len(dummy)]
            del path_l[0:i]
        if k <= i:
            dummy = path_l[i:j]
            for m in range(len(dummy)):
                path_l.insert(k+m, dummy[m])
            del path_l[i+len(dummy):j+len(dummy)]
        if k >= j:
            dummy = path_l[i:j]
            for m in range(len(dummy)):
                path_l.insert(k+m,dummy[m])
            del path_l[i:j] 
            
        path[d] = np.array(path_l)
        
    return path
 
       
#--------------------------------------------------------------------------|
 
       
def anneal_BRev_distance(Tem, path, nstep):
    '''
    Executes the annealing procedure by using block reverse moves,
    and the minimization of the distance as the move acceptance criterion.
    Parameters:
        Tem: temperatures' array.
        path: coordinates of the cities.
        nstep: number of moves for a single iteration.
    Returns:
        path: updated path.
        accelist: list containing the acceptance rate for each temperature.   
    '''
    accelist = []
    old=length(path)
    N = len(path[0])
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
 
       
def anneal_BRev_Metropolis(Tem, path, nstep):
    '''
    Executes the annealing procedure by using block reverse moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        Tem: temperatures' array.
        path: coordinates of the cities.
        nstep: number of moves for a single iteration.
    Returns:
        path: updated path.
        accelist: list containing the acceptance rate for each temperature.   
    '''
    accelist = []
    old=length(path)
    N = len(path[0])
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
            
            if np.exp(-(new-old)/Tem[i]) > rand.rand():
                old = new
                acce += 1
                new = length(path)
            else:
                path = breverse(r1, r2, path)
                
        accelist.append(100*float(acce)/nstep)
        
    return path, accelist

        
def anneal_swap_Metropolis(Tem, path, nstep):
    '''
    Executes the annealing procedure by using swap moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        Tem: temperatures' array.
        path: coordinates of the cities.
        nstep: number of moves for a single iteration.
    Returns:
        path: updated path.
        accelist: list containing the acceptance rate for each temperature.   
    '''
    accelist = []
    old = length(path)
    N =len(path[0])
    for i in range(len(Tem)):
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            r1 = int(N*rand.rand())
            r2 = int(N*rand.rand())
            while r1 > (N-1) and r2 > (N-1) : # In order to avoid IndexError
                r1 = int(N*rand.rand())
                r2 = int(N*rand.rand())
            path = swap(r1, r2, path)
            new = length(path)
            
            if np.exp(-(new-old)/Tem[i]) > rand.rand():
                old = new
                acce += 1
                new = length(path)
            else:
                path = swap(r2, r1, path)   
                
        accelist.append(100*float(acce)/nstep)
        
    return path, accelist


def anneal_swap_distance(Tem, path, nstep):
    '''
    Executes the annealing procedure by using swap moves,
    and the minimization of the distance as the move acceptance criterion.
    Parameters:
        Tem: temperatures' array.
        path: coordinates of the cities.
        nstep: number of moves for a single iteration.
    Returns:
        path: updated path.
        accelist: list containing the acceptance rate for each temperature.   
    '''
    accelist = []
    old = length(path)
    N = len(path[0])
    for i in range(len(Tem)):
        acce=0; ii=0
        while ii <= nstep: 
            ii += 1
            r1 = int(N*rand.rand())
            r2 = int(N*rand.rand())
            while r1 > (N-1) and r2 > (N-1) : # In order not to get and IndexError
                r1 = int(N*rand.rand())
                r2 = int(N*rand.rand())
            path = swap(r1, r2, path)
            new = length(path)
            
            if np.exp(-(new-old)/Tem[i]) >= 1:
                old = new
                acce += 1
                new = length(path)
            else:
                path = swap(r2, r1, path)
        
        accelist.append(100*float(acce)/nstep)
        
    return path, accelist


def anneal_PG_Metropolis(Tem, path, nstep):
    '''
    Executes the annealing procedure by using prune and graft moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        Tem: temperatures' array.
        path: coordinates of the cities.
        nstep: number of moves for a single iteration.
    Returns:
        path: updated path.
        accelist: list containing the acceptance rate for each temperature.   
    '''
    accelist = []
    old = length(path)
    N = len(path[0])
    for i in range(len(Tem)):
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            r1 = int(N*rand.rand())
            r2 = int(N*rand.rand())
            r3 = int(N*rand.rand())
            while r1 > (N-1) and r2 > (N-1) and r3 > (N-1) : # In order to avoid IndexError
                r1 = int(N*rand.rand())
                r2 = int(N*rand.rand())
                r3 = int(N*rand.rand())
            
            backup = path        
            path = prunegraft(r1, r2, r3, path)
            new = length(path)
            
            if np.exp(-(new-old)/Tem[i]) > rand.rand():
                old = new
                acce += 1
                new = length(path)
            else:
                path = backup
        
        accelist.append(100*float(acce)/nstep)
        
    return path, accelist


def anneal_PG_distance(Tem, path, nstep):
    '''
    Executes the annealing procedure by using prune and graft moves,
    and the Metropolis algorithm as the move acceptance criterion.
    Parameters:
        Tem: temperatures' array.
        path: coordinates of the cities.
        nstep: number of moves for a single iteration.
    Returns:
        path: updated path.
        accelist: list containing the acceptance rate for each temperature.   
    '''
    accelist = []
    old = length(path)
    N = len(path[0])
    for i in range(len(Tem)):
        acce = 0; ii = 0
        while ii <= nstep:
            ii += 1
            r1 = int(N*rand.rand())
            r2 = int(N*rand.rand())
            r3 = int(N*rand.rand())
            while r1 > (N-1) and r2 > (N-1) and r3 > (N-1) : # In order to avoid IndexError
                r1 = int(N*rand.rand())
                r2 = int(N*rand.rand())
                r3 = int(N*rand.rand())
            
            backup = path
            path = prunegraft(r1, r2, r3, path)
            new = length(path)
            
            if np.exp(-(new-old)/Tem[i]) >= 1:
                old = new
                acce += 1
                new = length(path)
            else:
                path = backup
            
        accelist.append(100*float(acce)/nstep)
        
    return path, accelist
