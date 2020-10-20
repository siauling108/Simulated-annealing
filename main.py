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

file1 = config.get('files','distances')
file2 = config.get('files','path')
file3 = config.get('files','Tem')
file4 = config.get('files','tot_accept')

while True:
    choice = input('''Choose the move type + the acceptance criterion to be used in
                   the annealing procedure:\n
                       BD = block reverse, distance;\n
                       BM = block reverse, Metropolis;\n
                       SM = swap ,Metropolis;\n
                       SD = swap, distance;\n
                       PM = prune and graft, Metropolis;\n
                       PD = prune and graft, distance.\n''')
    if choice in ['BD', 'BM', 'SM', 'SD', 'PM', 'PD']:
        break

path = functions.travel(N)
start=functions.length(path)
print('Initial distance:', round(start, 4))
distances=[start]
Tem = functions.Temp(T, T_min, alpha) #I need it in order to def. tot_acceptance
tot_acceptance = np.zeros((iterations, len(Tem)))

for ii in range(iterations):
    if choice == 'BD': 
        path, accelist=functions.anneal_BRev_distance(N, Tem, path, nstep)
    if choice == 'BM':
        path, accelist=functions.anneal_BRev_Metropolis(N, Tem, path, nstep) 
    if choice == 'SM':
        path, accelist=functions.anneal_swap_Metropolis(N, Tem, path, nstep)
    if choice == 'SD':
        path, accelist=functions.anneal_swap_distance(N, Tem, path, nstep) 
    if choice == 'PM':
        path, accelist=functions.anneal_PG_Metropolis(N, Tem, path, nstep)
    if choice == 'PD':
        path, accelist=functions.anneal_PG_distance(N, Tem, path, nstep)
    
    tot_acceptance[ii][:]=accelist
    distances.append(functions.length(path))
    print('Current distance:', round(distances[ii+1], 4),', Iteration:',ii+1)

np.save(file1, distances)
np.save(file2, path)
np.save(file3, Tem)
np.save(file4, tot_acceptance)
    
print('Final distance:',round(distances[-1], 4))
print("Percentage decrease of the distance", round((distances[0]-distances[-1])/distances[0]*100, 4), "%")
