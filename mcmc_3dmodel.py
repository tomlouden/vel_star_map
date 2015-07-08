# -*- coding: utf-8 -*-
from model3d import *
import numpy as np
import os
import pickle
from sampling_tools import sample_with_replacement, sample_without_replacement
import emcee
from random import randint

def main():

  bootstrap = False
#  x = [2.26339698e+00,-2.05583762e+00,8.55614089e+00,-7.98945634e-02,1.07223682e+00,9.48587317e+44,3.66446647e+00]
#  best_fit={'atm_radius':0.15}

  #x = [2.26339698e+00,-2.05583762e+00,8.55614089e+00,-7.98945634e-02,0.275,1.07223682e+00,9.48587317e+44,3.66446647e+00]
  #x = [1.66339698e+00,-0.55583762e+00,8.55614089e+00,-7.98945634e-02,0.15,1.07223682e+00,9.48587317e+44,3.66446647e+00]

#  for when limb darkening is fixed
  x = [7.66635193,-0.09572382,0.21485677,0.58513335,1.21818912e+45,3.66011491]
  best_fit={'ld1':0.51374991,'ld2':0.18869174}

  write_file = '250w_fixedld_synth_spectrum.txt'
  resume = False

  nwalkers = 250
  ndim = len(x)
  nproc = 4

  mu = 1.0
  sigma = 0.1


  p0 = []

  if resume == True:
    data = np.loadtxt(write_file)
    p0 = data[-20*16:]
    seed = [randint(0,20*16) for p in range(0,nwalkers)]
  else:
    for i in range(0,nwalkers):
      lp = -np.inf    
      while lp == -np.inf:
        p = np.random.normal(mu, sigma, ndim)
        lp = lnprior(p*x)
      print lp
      p0 += [p*x]
    p0 = array(p0)

  sampler = emcee.EnsembleSampler(nwalkers, ndim, initialise_walker, args=[nproc,bootstrap,best_fit])

#  f = open(write_file, "w")
#  f.close()

  for result in sampler.sample(p0, iterations=10000, storechain=False):
    position = result[0]
    f = open(write_file, "a")
    for k in range(position.shape[0]):
      print position[k]
#      f.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
      outwrite = ''
      for i in range(0,len(position[k])):
        outwrite += str(position[k][i])+' ' 
      f.write(outwrite+'\n')

    f.close()

  quit()

  initialise_walker(x,bootstrap=bootstrap,best_fit=best_fit)

def initialise_walker(x,nproc,bootstrap=False,best_fit=False):

  sodium_d = [5889.950,5895.924]
  midtransit = 54341.056842
  period =  2.21857312
  planet_K = 154
  star_K = 0.20196
  line_centers = sodium_d

  if bootstrap == True:
    index, index2 = sample_with_replacement(20,20)
    index = index + 9
    index = np.sort(index)
  else:
    index = np.arange(9,29)

  plotting = False

#  lp = lnprior(x)

  lp = 0.0

  if lp == -np.inf:
    print 'fail!'
    print x[0],x[1]
    print lp
    return lp  

  lnprob = calc_profile(x,index,midtransit,period,planet_K,star_K,line_centers,plotting=plotting,nproc=nproc,best_fit=best_fit)

  lnprob = lnprob + lp

  return lnprob

def lnprior(x):
  limb = x[0] + x[1]
  if limb < 1:
    return -np.inf
  else:
    return 0


def calc_profile(x,clipping_index,midtransit,period,planet_K,star_K,line_centers,plotting=True,nproc=4,best_fit=False):

  # load the master spectrum
  master_dat = np.loadtxt('sodium_spectrum.dat')
  master_dat = np.loadtxt('synthetic_spectrum.dat')

  master_wvl = master_dat[:,0]
  master_flux = master_dat[:,1]

# the old way - probably contaminated by my mistake with the rvs

  #dir_old = 'sodium_diff_bin/'

  dir_old = 'no_poly_correct_rv/'
  dir = 'no_poly_correct_rv/'
  dir_old = 'everything/'
  dir_old = 'errors/'

  dir = 'model_rv/'
  dir = 'test_rv/'
  dir = 'everything/'
  #dir = 'errors/'
  #dir = 'ccf_rvs/'
  #dir = 'diff_only/'


  files = os.listdir(dir)


  sodium_d = [5889.950,5895.924]

  time = []
  wvl = []
  spectra = []
  spectra_errors = []
  spectra2 = []
  spectra_errors2 = []

  min_wvl = 5887
  max_wvl = 5899

  min_wvl = 5869
  max_wvl = 5909

  telluric = pickle.load(open('example_telluric.p','rb'))

  for f in files:
    time += [float(f[:-3])]
    data = pickle.load(open(dir+f,'rb'))
    data2 = pickle.load(open(dir_old+f,'rb'))
    index = [(data['wvl'] > min_wvl) & (data['wvl'] < max_wvl)]

    wvl += [data['wvl'][index]]
    spectra += [data['spec'][index]]
    spectra_errors += [data['spec_error'][index]]
    spectra2 += [data2['spec'][index]]
    spectra_errors2 += [data2['spec_error'][index]]

  min_wvl = 5887
  max_wvl = 5899

  min_wvl = wvl[0][np.argmin(wvl[0])+2]
  max_wvl = wvl[0][np.argmax(wvl[0])-2]

  min_wvl = 5887 - 9
  max_wvl = 5899 + 9

  time = array(time)

  spectra = array(spectra)
  spectra_errors = array(spectra_errors)
  wvl = array(wvl)

  # don't try to cross correlate... is nasty 
  #rv = measure_rvs(wvl,spectra,min_vel=-30000,max_vel=30000,steps=10)
  #plot(time,rv,'bo')
  #show()

  index = [np.argsort(time)]
  #index = [arange(0,len(time))]
  time = time[index]
  spectra = spectra[index]
  spectra_errors = spectra_errors[index]
  #spectra_errors = spectra_errors[index]*0 + 1.0
  wvl = array(wvl[index])

  wvl = array(wvl[0])

  #stellar_rv = -201.96*sin(2*pi*(midtransit-time)/period)
  #spectra = align_spectra(wvl,spectra,stellar_rv)
  #spectra_errors = align_spectra(wvl,spectra_errors,stellar_rv)

  spectra = spectra[clipping_index]
  spectra_errors = spectra_errors[clipping_index]
  time = time[clipping_index]

  index = [np.argsort(time)]
  #index = [arange(0,len(time))]
  time = time[index]
  spectra = spectra[index]
  spectra_errors = spectra_errors[index]

  c = 299792458.0

  velocity_line1 = -((wvl - sodium_d[0])/sodium_d[0])*c
  velocity_line2 = -((wvl - sodium_d[1])/sodium_d[1])*c

  combined_vel = array(list(velocity_line1) + list(velocity_line2))
  vel_index = np.argsort(combined_vel)
  combined_vel = combined_vel[vel_index]

#  show()
  
  resamp_factor = 10
  dummy_wvl = np.linspace(min(wvl),max(wvl),len(wvl)/resamp_factor)


  # time to try the 3d model

  prior_3d = []

  for i in range(0,len(spectra)):
    spectra[i][0] = 1.0
    spectra[i][-1] = 1.0

  gamma_0 = 0.94426940
  gamma_1 = -0.41811691
  planet_K = 154
  vel_offset = 2.3
  atm_radius = 0.1
  fwhm = 0.5
  ratio = 300

  diff = fit_model_3d(x,time,period,planet_K,star_K,midtransit,spectra,spectra_errors,wvl,line_centers,plotting=False,nproc=nproc,best_fit=best_fit)

  elements = len(diff)

  chi2 = sum(diff**2)/elements

  print chi2

  elements = len(spectra)*len(spectra[0])

  sn = spectra_errors.reshape((elements))

  lnprob = -0.5*sum((diff**2/sn**2) + np.log(2.0*pi*sn**2))
 
  print lnprob

  return lnprob

main()
