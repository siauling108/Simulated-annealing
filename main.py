# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 09:53:20 2020

@author: marco
"""
import numpy as np
import configparser
#import sys
#from sys import argv
import functions


config = configparser.ConfigParser()
config.read('input.txt')
#config.read(sys.argv[1])

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
    choice = input('''Choose the method + the acceptance criterion to be used in
                   the annealing procedure:\n
                       BD = block reverse, distance;\n
                       BM = block reverse, Metropolis;\n
                       SM = swap ,Metropolis;\n
                       SD = swap, distance;\n
                       PM = prune and graft, Metropolis;\n
                       PD = prune and graft, distance.\n''')
    if choice in ['BD', 'BM', 'SM', 'SD', 'PM', 'PD']:
        break

city = functions.travel(N)
citynum = list(range(N))
start=functions.length(N, city, citynum)
print('Initial distance:',start)
distances=[start]
T_len = functions.Tem_length(T, T_min, alpha) #I need it in order to def. tot_acceptance
tot_acceptance = np.zeros((iterations, T_len))

for ii in range(iterations):
    if choice == 'BD': 
        citynum, accelist, Tem=functions.anneal_BRev_distance(N, alpha, T, city, T_min, nstep, citynum)
    if choice == 'BM':
        citynum, accelist, Tem=functions.anneal_BRev_Metropolis(N, alpha, T, city, T_min, nstep, citynum) 
    if choice == 'SM':
        citynum, accelist, Tem=functions.anneal_swap_Metropolis(N, alpha, T, city, T_min, nstep, citynum) #BOH
    if choice == 'SD':
        citynum, accelist, Tem=functions.anneal_swap_distance(N, alpha, T, city, T_min, nstep, citynum) 
    if choice == 'PM':
        citynum, accelist, Tem=functions.anneal_PG_Metropolis(N, alpha, T, city, T_min, nstep, citynum)
    if choice == 'PD':
        citynum, accelist, Tem=functions.anneal_PG_distance(N, alpha, T, city, T_min, nstep, citynum)

    tot_acceptance[ii][:]=accelist
    distances.append(functions.length(N, city, citynum))
    print('Current distance:',distances[ii+1],', Iteration:',ii+1)

path = functions.get_path(N, city, citynum)

np.save(file1, distances)
np.save(file2, path)
np.save(file3, Tem)
np.save(file4, tot_acceptance)
    
print('Final distance:',distances[-1])
print("Percentage decrease of the distance", (distances[0]-distances[-1])/distances[0]*100., "%")
