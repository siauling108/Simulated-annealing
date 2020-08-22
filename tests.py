# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 09:57:11 2020

@author: marco
"""
import hypothesis
from hypothesis import strategies as st
from hypothesis import given, assume, settings
import hypothesis.extra.numpy as xnp
import numpy as np
#import numpy.random as rand
import functions
import configparser


config = configparser.ConfigParser()
config.read('inputs.txt')

N = config.getint('parameters', 'N') 
iterations = config.getint('parameters', 'iterations')
nstep = config.getint('parameters', 'nstep')
T_min = config.getfloat('parameters', 'T_min')
T = config.getfloat('parameters', 'T')
alpha = config.getfloat('parameters', 'alpha')


@given(N=st.integers(1, N))
@settings(max_examples = 1)
def test_city(N):
    city = functions.travel(N)
    assert len(city) == 2


@given(M=st.integers(2,N), city=xnp.arrays(np.float,(2,N),st.floats(0.01,1)), citynum=st.lists(st.integers(),min_size=2,max_size=N))
@settings(max_examples = 1)
def test_leng_positive(M, city, citynum):    #doesn't work for now
    #assume(len(citynum)>1)
    #assume(all(i > 0 for i in city))        
    distance = functions.length(M, city, citynum)
    assert distance > 0
