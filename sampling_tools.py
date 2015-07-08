# -*- coding: utf-8 -*-
from numpy.random import randint
import numpy as np

def sample_with_replacement(n,s1):
  l2 = []
  l1 = []

  sample = list(np.arange(0,n))

  for i in range(0,s1):
    r = randint(0,len(sample))
    l1 += [sample[r]]

  for i in range(0,n-s1):
    r = randint(0,len(sample))
    l2 += [sample[r]]

  l1 = np.array(l1)
  l2 = np.array(l2)
  
  return l1, l2

def sample_without_replacement(n,s1):
  l2 = list(np.arange(0,n))
  l1 = []
  for i in range(0,s1):
    r = randint(0,len(l2))
    l1 += [l2[r]]
    l2.remove(l2[r])
  
  l1 = np.array(l1)
  l2 = np.array(l2)
  return l1, l2