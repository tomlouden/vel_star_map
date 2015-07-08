# -*- coding: utf-8 -*-
from pysynphot import observation
from pysynphot import spectrum
from numpy import *
import numpy as np
#from pylab import *
def measure_rvs(wvl,spectra,min_vel=-300,max_vel=300,steps=100):

  # make a master spectrum
  master_spectra = np.sum(spectra, axis = 0)
  master_spectra = master_spectra/np.median(master_spectra)

  rvs = []
  for i in range(0,len(spectra)):
    rv = find_rv(wvl,spectra[i]/np.median(spectra[i]),wvl,master_spectra,min_vel=min_vel,max_vel=max_vel,steps=steps)
    rvs += [rv]
    print rv
  rvs = np.array(rvs)
  return rvs

def rebin_spec(wave,specin,wavnew):
  spec = spectrum.ArraySourceSpectrum(wave=wave , flux=specin)
  f = np.ones(len(wave))
  filt = spectrum.ArraySpectralElement(wave,f,waveunits='angstrom')
  obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
  return obs.binflux

def super_sample(in_array,factor):

  super_sampled = []
  for i in in_array:
    super_sampled += [i]*factor

  super_sampled = array(super_sampled)
  return super_sampled

def redshift(wave,specin,velocity):
  
  c = 299792458.0

  #appy redshift
  wave_shift = wave*(1.0 - (velocity/c))

  # regrid the spectrum onto the original wavelength scale
  specout = rebin_spec(wave_shift,specin,wave)

  return specout

def find_rv(wvl,spec,ref_wvl,ref_spec,min_vel=-300,max_vel=300,steps=1000):
  import scipy.optimize as opt
  from scipy.interpolate import interp1d
  from pylab import *

  # make velocity grid
  vels = np.linspace(min_vel,max_vel,steps)

  # find the ccf
  ccf = calc_ccf(vels,wvl,spec,ref_wvl,ref_spec)

  #plot(vels,ccf)
  #show()

  # fit for minimum of the ccf
  #gauss_prior = array([1000,vels[argmin(ccf)],10000,1000])
  #final, cov = opt.leastsq(fit_gauss,gauss_prior,args=(vels,ccf))
  #model = gauss_model(final,vels)
  #plot(vels,ccf,'ro')
  #plot(vels,model,'b-')
  #show()
  #rv = final[1]

  # or do a cubic spline instead?
  #min_vel = vels[argmin(ccf)-10:argmin(ccf)+10]
  #min_ccf = ccf[argmin(ccf)-10:argmin(ccf)+10]
  #new_vel = np.linspace(min(min_vel),max(min_vel),1000)
  #f = interp1d(min_vel,min_ccf, kind='cubic')
  #rv = new_vel[argmin(f(new_vel))]

  # or, do something a bit cleverer, and actually fit for the minimum directly...
  final, cov = opt.leastsq(calc_ccf,vels[argmin(ccf)],args=(wvl,spec,ref_wvl,ref_spec))
  rv = final[0]

  return rv

def align_spectra(wvl,spectra,rv):
  for i in range(0,len(spectra)):
    spectra[i] = redshift(wvl,spectra[i],rv[i])
  return spectra

def fit_gauss(x,coord,data):
  model = gauss_model(x,coord)
  diff = model - data
  return diff

def gauss_model(x,coord):

  amp = x[0]
  line_cen = x[1]
  fwhm = abs(x[2])
  baseline = x[3]

  sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))

  gauss = np.exp(-0.5*((coord-line_cen)/sigma)**2.0)

  line = amp*gauss/max(gauss)

  model = baseline - line

  return model

def calc_ccf(vels,wvl,spec,ref_wvl,ref_spec):
  from pylab import *

  # rebin reference spectrum to observed spectrum
  ref_spec = rebin_spec(ref_wvl,ref_spec,wvl)

  ccf = []

  clipsize = 1

  for vel in vels:
    shift_spec = redshift(wvl,spec,vel)
    difference = (shift_spec[clipsize:-(clipsize+1)] - ref_spec[clipsize:-(clipsize+1)])/(ref_spec[clipsize:-(clipsize+1 )])
    ccf += [sum(difference**2)]
  ccf = np.array(ccf)

  #plot(wvl[10:-10],difference)
  #show()
  return ccf
