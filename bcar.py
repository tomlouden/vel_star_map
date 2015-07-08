# -*- coding: utf-8 -*-
from pylab import *
from scipy.stats import norm
from scipy.stats import uniform

def main():
  data = fakeit(1000,0)

  data = absorb_line(data,10)

  folded = boxcar(data,10)
  plot(data)
  plot(folded)
  show()

def absorb_line(data,fwhm):

  a = 1.0
  b = len(data)/2
  c = fwhm/2.35482
  x = np.arange(0,len(data))
  print a,b,c
  data -= a*exp(-((x-b)**2)/(2*c**2))
  return data


def boxcar(xaxis,data,n):
  output = array([mean(data[x:x+n]) for x in xrange(0,len(data)-n)])
  return xaxis[n/2:-n/2],output

def boxcar_detection(xaxis,data,n):
  output = array([data[x:x+n] for x in xrange(n,len(data)-n)])
  return xaxis[n:-n],output

def fakeit(n,z):
  length = array(np.arange(0,n))
  fake_curve = norm.rvs(size=n)
  for i in range(0,z):
    fake_curve += sin(length/uniform.rvs(0,n) + uniform.rvs(-2*pi,2*pi))
  return fake_curve
