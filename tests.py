# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 09:57:11 2020

@author: marco
"""
import hypothesis
from hypothesis import strategies as st
from hypothesis import given, assume
import functions
import configparser

config = configparser.ConfigParser()
config.read('input.txt')

N = config.getint('parameters', 'N') 
iterations = config.getint('parameters', 'iterations')
nstep = config.getint('parameters', 'nstep')
T_min = config.getfloat('parameters', 'T_min')
T = config.getfloat('parameters', 'T')
alpha = config.getfloat('parameters', 'alpha')


@given(N=st.integers(1, N), city=st.lists(), citynum=st.lists(min_size=0 , N))
def test_leng_positive(N, city, citynum):               #doesn't work for now
    distance = functions.length(N, city, citynum)
    assert distance > 0
