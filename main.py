# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 09:53:20 2020

@author: marco
"""
import numpy as np
import configparser
import functions


config = configparser.ConfigParser()
config.read('inputs.txt')

N = config.getint('parameters', 'N') 
iterations = config.getint('parameters', 'iterations')
nstep = config.getint('parameters', 'nstep')
T_min = config.getfloat('parameters', 'T_min')
T = config.getfloat('parameters', 'T')
alpha = config.getfloat('parameters', 'alpha')

method = config.get('parameters', 'method')
T_profile = config.get('parameters', 'T_profile')
travel = config.get('parameters', 'travel')

file1 = config.get('files','distances')
file2 = config.get('files','path')
file3 = config.get('files','Tem')
file4 = config.get('files','tot_accept')


#--------------------------------------------------------------|


if travel == 'travel_random':
    path = functions.travel_RND(N)

if T_profile == 'Temp_decr':
    Tem = functions.Temp_decr(T, T_min, alpha) 

start=functions.length(path)
distances=[start]
tot_acceptance = np.zeros((iterations, len(Tem)))

print('Method:', method)
print('Travel:', travel)
print('Initial distance:', round(start, 4))

for ii in range(iterations):
    if method == 'BD': 
        path, accelist=functions.anneal_BRev_distance(Tem, path, nstep)
    if method == 'BM':
        path, accelist=functions.anneal_BRev_Metropolis(Tem, path, nstep) 
    if method == 'SM':
        path, accelist=functions.anneal_swap_Metropolis(Tem, path, nstep)
    if method == 'SD':
        path, accelist=functions.anneal_swap_distance(Tem, path, nstep) 
    if method == 'PM':
        path, accelist=functions.anneal_PG_Metropolis(Tem, path, nstep)
    if method == 'PD':
        path, accelist=functions.anneal_PG_distance(Tem, path, nstep)
    
    tot_acceptance[ii][:]=accelist
    distances.append(functions.length(path))
    print('Current distance:', round(distances[ii+1], 4),', Iteration:',ii+1)


np.save(file1, distances)
np.save(file2, path)
np.save(file3, Tem)
np.save(file4, tot_acceptance)
    
print('Final distance:',round(distances[-1], 4))
print("Percentage decrease of the distance", round((distances[0]-distances[-1])/distances[0]*100, 4), "%")
