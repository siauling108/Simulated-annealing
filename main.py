# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 09:53:20 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import configparser
import sys
from sys import argv


config = configparser.ConfigParser()
config.read(sys.argv[1])

N = config.get('parameters', 'N')
iterations = config.get('parameters', 'iterations')
T_min = config.get('parameters', 'T_min')
T = config.get('parameters', 'T')

city = np.zeros((2,N))
lista = list(range(N))
#lista=list(lista)
Tem = []
accelist = []
for i in range(N):
    for j in (0,1):
        city[j][i] = rand.rand() #cities coordinates 
  