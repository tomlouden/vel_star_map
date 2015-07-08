# -*- coding: utf-8 -*-
import pickle
import os
from pylab import *
from numpy import *
import numpy as np
from spec_manip import *
#from correct_linearity import constrained_sodium_model, polynomial
import scipy.optimize as opt
from sampling_tools import sample_with_replacement, sample_without_replacement
from plot_spectrum import bin_curve
from matplotlib.ticker import MaxNLocator
#from mp_vel_prism import vel_prism
from mp_vel_prism_new import vel_prism
from pyspeckit.spectrum.models.inherited_voigtfitter import voigt
from model3d import *

size = 16
rcParams['font.family'] ='serif'
rcParams['font.serif'] ='Computer Modern Roman'
rcParams['text.usetex'] =True
rcParams['axes.labelsize'] =size
rcParams['xtick.labelsize'] =size
rcParams['ytick.labelsize'] =size
rcParams['legend.fontsize'] =size
rcParams['text.usetex'] =True
rcParams['figure.figsize'] =14.6,8.4

def main():

  nproc = 20

  sodium_d = [5889.950,5895.924]
  midtransit = 54341.056842
  period =  2.21857312
  planet_K = 154
  star_K = 0.20196
  line_centers = sodium_d

  index = arange(9,29)

  plotting = True
  
  #output = calc_profile(sort(index),midtransit,period,planet_K,star_K,line_centers,plotting=plotting,nproc=1)
  #quit()

  for x in range(0,10000):
#    index, index2 = sample_without_replacement(20,15)
    index, index2 = sample_with_replacement(20,20)
    index = index + 9
    output = calc_profile(index,midtransit,period,planet_K,star_K,line_centers,plotting=False,nproc=nproc)
    with open('bootstrap_ldfree.txt','a') as outfile:
      outwrite = ''
      for i in range(0,len(output)):
	outwrite += str(output[i])+' ' 
      outfile.write(outwrite+'\n')


def calc_profile(clipping_index,midtransit,period,planet_K,star_K,line_centers,plotting=True,nproc=4):

  best_fit = [-1.14601379e+02,3.67424540e+00,6.97253518e-03,4.89808018e-01]

  # load the master spectrum
  master_dat = np.loadtxt('sodium_spectrum.dat')
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
    bin_wvl, bin_combined, bin_combined_error = bin_curve(wvl[-1],spectra[-1],spectra_errors[-1],1)
    bin_wvl2, bin_combined2, bin_combined_error2 = bin_curve(wvl[-1],spectra2[-1],spectra_errors2[-1],1)

    #p0 = 0.1572
    #radiustotal = (1.0 + p0)/8.92
    #gamma_0 = 0.94426940
    #gamma_1 = -0.41811691
    #inc = 85.68
    #star_vsini = 3.1
    #planet_K = 154
    #vel_offset = 2.3
    #atm_radius = 0.1
    #input = [p0,radiustotal,gamma_0,gamma_1,inc,star_vsini,planet_K,vel_offset,atm_radius]
    #fwhm = 0.5
    #ratio = 300
    #data_x = array([(time[-1]-midtransit)/period])
    #master_dat = loadtxt('sodium_spectrum.dat')
    #model_wvl = array(master_dat[:,0])
    #master_flux = array(master_dat[:,1])
    #profile = {'wvl':model_wvl,'spectrum':master_flux}
    #results = vel_prism(data_x, input, profile, time, spot_data=False, type_data='FLUX',plotting=False,result_wvl=wvl[-1])
    #errorbar(bin_wvl,bin_combined,yerr=bin_combined_error,fmt='r.')
    #plot(wvl[-1],results[0])
    #show()

    #plot(bin_wvl,bin_combined,'.')
    #errorbar(bin_wvl,bin_combined2,yerr=bin_combined_error2,fmt='r.')
    #errorbar(bin_wvl,bin_combined,yerr=bin_combined_error,fmt='r.')
    #errorbar(bin_wvl,bin_combined2,yerr=bin_combined_error2,fmt='g.')
    #plot(bin_wvl,1.0+bin_combined2-bin_combined,'r')

    #plot(telluric['wvl'],1.0 + (telluric['spectra']-1.0)/39,'k')
    #plot(master_wvl,1+((master_flux/median(master_flux))-1.0)/40,'b')

  #show()

  #plot(wvl[0],sum(spectra,axis=0),'bo')
  #show()

  min_wvl = 5887
  max_wvl = 5899

  min_wvl = wvl[0][argmin(wvl[0])+2]
  max_wvl = wvl[0][argmax(wvl[0])-2]

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

  index = [argsort(time)]
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

  index = [argsort(time)]
  #index = [arange(0,len(time))]
  time = time[index]
  spectra = spectra[index]
  spectra_errors = spectra_errors[index]
  #spectra_errors = spectra_errors[index]*0 + 1.0

  #print abs(time - 54341.0568358)*24 > 0.48613901131

  # sanity check to make sure we're alligning properly
  #vel_model = 0*best_fit[1] + -1*best_fit[0]*sin(2.0*pi*(time-midtransit)/period)
  #t = (time - midtransit)*24
  #norm = max(abs(t))
  #rv = 1000*0.0*sin((pi*t/norm)/2.0)
  #rv = 1000*0.0*sin((pi*((t/norm)+0.15))/2.0)

  #vel_model = 1000*vel_model + rv
  #spectra = align_spectra(wvl,spectra.copy(),vel_model)

  c = 299792458.0

  velocity_line1 = -((wvl - sodium_d[0])/sodium_d[0])*c
  velocity_line2 = -((wvl - sodium_d[1])/sodium_d[1])*c

  combined_vel = array(list(velocity_line1) + list(velocity_line2))
  vel_index = argsort(combined_vel)
  combined_vel = combined_vel[vel_index]

#  show()
  
  resamp_factor = 10
  dummy_wvl = linspace(min(wvl),max(wvl),len(wvl)/resamp_factor)


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

  prior_3d = [gamma_0,gamma_1,planet_K,vel_offset,atm_radius,fwhm,ratio]

  # for when you're fixed to 0.1
  prior_3d = [1.65308844e+00,-1.04812914e+00,1.00491621e+02,5.71777342e+00,1.00000000e-01,2.13167534e+00,-4.56788205e+07]

  # for when you're fixed to 0.2
  prior_3d = [1.97972829e+00,-1.76607321e+00,1.19430590e+02,4.84616011e+00,2.00000000e-01,1.20117236e+00,-2.75251623e+14]

  # added the east/west offset split
  prior_3d = [1.98496666e+00,-1.78689117e+00,1.54000000e+02,7.35833821e+00,1.49665948e+00,3.00000000e-01,7.71559205e-01,3.68166783e+28]


  # careful! from now on we allow K to be fixed. (could make it a gaussian prior?)
  prior_3d = [1.98496666e+00,-1.78689117e+00,7.35833821e+00,1.49665948e+00,3.00000000e-01,7.71559205e-01,3.68166783e+28]
  best_fit=False

  # careful! from now on we leave radius fixed until we figure out the genetic algorithm
  prior_3d = [1.98496666e+00,-1.78689117e+00,7.35833821e+00,1.49665948e+00,7.71559205e-01,3.68166783e+28]
  prior_3d = [1.95876320e+00,-1.73341318e+00,6.56348290e+00,2.29320705e+00,9.45055108e-01,-6.48816563e+33]

  prior_3d = [2.73989049e+00,-2.86942945e+00,6.74199449e+00,1.99626115e+00,9.52902509e-01,9.04716892e+40]
  best_fit={'atm_radius':0.275}

  #prior_3d = [2.23924678e+00,-2.02719249e+00,5.48955581e+00,3.45310090e+00,2.39205669e+00,6.38920070e+39]
  #best_fit={'atm_radius':0.1}

  # now we've added saturation! line profiles much more realistic

  prior_3d = [2.23924678e+00,-2.02719249e+00,5.48955581e+00,3.45310090e+00,2.39205669e+00,6.38920070e+39,2.0]
  prior_3d = [2.26339698e+00,-2.05583762e+00,8.55614089e+00,-7.98945634e-02,1.07223682e+00,9.48587317e+44,3.66446647e+00]
  prior_3d = [1.8,-0.8,8.55614089e+00,-7.98945634e-02,1.07223682e+00,9.48587317e+44,3.66446647e+00]
  best_fit={'atm_radius':0.15}

  prior_3d = [0.51374991, 0.18869174, 7.3171305975049998, -0.083571968753300002, 0.44724670399400002, 1.24645711824e+45, 4.2957133770149998]
  best_fit={'atm_radius':0.235380980973}

#  prior_3d = [7.3171305975049998, -0.083571968753300002, 0.44724670399400002, 1.24645711824e+45, 4.2957133770149998]
#  best_fit={'ld1':0.51374991,'ld2':0.18869174,'atm_radius':0.235380980973}

  output = opt.leastsq(fit_model_3d,prior_3d,args=(time,period,planet_K,star_K,midtransit,spectra,spectra_errors,wvl,line_centers,plotting,nproc,best_fit),full_output=True)

  final = output[0]
  cov = output[1]

  return final

  fit_model_3d(final,time,period,planet_K,star_K,midtransit,spectra,spectra_errors,wvl,line_centers,True,1,best_fit=best_fit)

  #linetype = 'lorentz'
  #linetype = 'gauss'
  linetype = 'voigt'

  if linetype == 'voigt':
    prior = [130,5,5.0e-3,0.5,0.5]
  else:
    prior = [13,0.5,5.0e-3,0.5]

  #poly_prior = [1,1e-5,1e-9,1e-13,1e-16]
  poly_prior = [1]
  for i in range(0,len(spectra)):
    prior += poly_prior
  len_pol = len(poly_prior)

  #output = opt.leastsq(sin_func,prior,args=(time,period,midtransit,spectra,spectra_errors,wvl,len_pol,False,linetype),full_output=True)

  free_prior = [0.5,5.0e-3,300] + [10]*20
  poly_prior = [1]
  for i in range(0,len(spectra)):
    free_prior += poly_prior
  len_pol = len(poly_prior)

  output = opt.leastsq(free_func,free_prior,args=(time,period,midtransit,spectra,spectra_errors,wvl,len_pol,plotting,linetype),full_output=True)

  #output = opt.minimize(sin_func,prior,args=(time,period,midtransit,spectra,spectra_errors,wvl,len_pol))

  final = output[0]
  cov = output[1]

  print final

  vel_model = final[3:23]

  print 'mean, median, std'
  print 'ingress ',mean(vel_model[0:5]),median(vel_model[0:5]),std(vel_model[0:5])
  print 'full ',mean(vel_model[5:15]),median(vel_model[5:15]),std(vel_model[5:15])
  print 'egress ',mean(vel_model[15:20]),median(vel_model[15:20]),std(vel_model[15:20])

  plot(vel_model,'bo')

  show()

  quit()

  # let's see what we have

  if plotting != True:
    return final

  print final

  model = sin_model(final,time,period,midtransit,wvl,len_pol,linetype)

  # right, we need to get in here and remove all the baselines
  # 

  baseline_vals = final.copy()

  baseline_vals[2] = 0 

  #uncomment this if you want to see the individual plots
  #sin_func(final,time,period,midtransit,spectra,spectra_errors,wvl,len_pol,plotting=True,linetype=linetype)

  model = array(model)
  baselines =   sin_model(baseline_vals,time,period,midtransit,wvl,len_pol,linetype,baseline=True)
  spectra = spectra/baselines

  npoints = len(model)*len(model[0])
  nfree = npoints - 1 - len(final)
  chi2 = sum(((spectra-model)/spectra_errors)**2.0)
  BIC = chi2 + len(final)*log(npoints)
  Rchi2 = chi2 / nfree
  try:
    s_errors = sqrt(abs(cov)*Rchi2)
    fit_errors = [s_errors[i][i] for i in range(0,len(s_errors))]
  except:
    fit_errors = final

  print (final[0]),'+-',fit_errors[0],'km/s projected orbital velocity'
  print (final[1]),'+-',fit_errors[1],'km/s offset velocity'
  print abs(final[2]),'+-',fit_errors[2],'Fractional depth'
  print abs(final[3]),'+-',fit_errors[3],'FWHM'

  vel_model2 = final[1] + final[0]*sin(2.0*pi*(midtransit-time)/period)
  vel_model2 = final[1] + 154*sin(2.0*pi*(midtransit-time)/period)

  print vel_model2[0], 'velocity of frame'
  print vel_model2[0]-final[1], 'orbital velocity of frame'

  vel_model = 0*final[1] + final[0]*sin(2.0*pi*(midtransit-time)/period)
  vel_model2 = final[1] + 0*final[0]*sin(2.0*pi*(midtransit-time)/period)

  vel_model = 0*final[1] + 154*sin(2.0*pi*(midtransit-time)/period)
  vel_model2 = final[1] + 0*final[0]*sin(2.0*pi*(midtransit-time)/period)

  d1 = 5895.924
  d2 = 5889.950
  velocity = 10e3
  c = 3e8
  alligned = align_spectra(wvl,spectra.copy(),vel_model)


  sod_model = [final[2],final[1],final[3],1.0]

  if linetype == 'voigt':
    sod_model = [final[2],final[1],final[3],final[4],1.0]

  smodel = constrained_sodium_model(sod_model,wvl,wvl.copy(),linetype)


  make_velspace_plot(alligned,smodel,wvl,d1,d2)

  fig, (ax1,ax2) = plt.subplots(2, sharex=True)

  ax1.minorticks_on()
  ax2.minorticks_on()

  t = (time - midtransit)*24
  norm = max(abs(t))
  rv = 1000*0.0*sin((pi*t/norm)/2.0)
  #rv = 1000*0.0*sin((pi*((t/norm)+0.15))/2.0)

  print (max(vel_model) - min(vel_model)), 'change from start to end of transit'

  vel_model2 = 1000*vel_model2
  vel_model = 1000*vel_model + rv

  sodium_d = [5889.950,5895.924]
  linewidth = 3.0


  print (max(vel_model) - min(vel_model))/1000, 'change from start to end of transit'

  print smodel

  ax1.set_ylim(0.0,1.2) 
  ax1.set_yticks([0.2,0.4,0.6,0.8,1.0])
  ax1.plot(master_wvl,master_flux/median(master_flux),'k-')
  narrow_region = [5887,5899]
  #ax1.axvline(x=narrow_region[0],color='k',ls='--')
  #ax1.axvline(x=narrow_region[1],color='k',ls='--')

  ax1.axvline(x=sodium_d[0]-linewidth/2.0,color='k',ls=':')
  ax1.axvline(x=sodium_d[0]+linewidth/2.0,color='k',ls=':')
  ax1.axvline(x=sodium_d[1]-linewidth/2.0,color='k',ls=':')
  ax1.axvline(x=sodium_d[1]+linewidth/2.0,color='k',ls=':')

  shifted_d1 = d1 - d1*velocity/c
  shifted_d2 = d2 - d2*velocity/c

  print shifted_d1
  print shifted_d2

  #ax1.axvline(x=shifted_d1,color='k',ls=':')
  #ax1.axvline(x=shifted_d2,color='k',ls=':')

  ax1.plot(telluric['wvl'],1.1 + 3*(telluric['spectra']-1.0)/39,'r')
  ax1.set_ylabel('Normalised intensity',labelpad=20)

  alligned = align_spectra(wvl,spectra.copy(),vel_model)

  #spectra = spectra/array([std(spectra,axis=1)]*2537).T

  #alligned = spectra.copy()

  #combined = median(alligned,axis=0)
  #combined = mean(alligned,axis=0)
  #errors = std(alligned,axis=0)/sqrt(len(spectra))
  #combined = rebin_spec(wvl,combined,dummy_wvl)
  #errors = rebin_spec(wvl,errors,dummy_wvl)/sqrt(resamp_factor)

  scale = sqrt(3.0)

  bin_wvl, bin_combined, bin_combined_error = bin_curve(wvl,median(alligned,axis=0),scale*std(alligned,axis=0)/sqrt(len(spectra)),resamp_factor)
  #bin_wvl, bin_combined, bin_combined_error = bin_curve(wvl,median(alligned,axis=0),std(alligned,axis=0)/sqrt(len(spectra)),resamp_factor)
  bin_wvl, bin_smodel, bin_smodel_error = bin_curve(wvl,smodel,std(alligned,axis=0)/sqrt(len(spectra)),resamp_factor)

  #ax2.errorbar(dummy_wvl[1:-2],combined[1:-2],yerr=errors[1:-2],fmt='bo')

  #ylim(0.992-1.0,1.003-1.0) 
  #ax2.plot(wvl,smodel-smodel,'r-')
  #ax2.errorbar(bin_wvl[1:-2],bin_combined[1:-2]-bin_smodel[1:-2],yerr=bin_combined_error[1:-2],fmt='k.')

  ax2.set_ylim(0.987,1.005) 
  ax2.set_ylim(0.987,1.010) 
  ax2.plot(wvl,smodel,'r-')
  ax2.errorbar(bin_wvl[5:-5],bin_combined[5:-5],yerr=bin_combined_error[5:-5],fmt='k.')

  #axvline(x=narrow_region[0],color='k',ls='--')
  #axvline(x=narrow_region[1],color='k',ls='--')

  axvline(x=sodium_d[0]-linewidth/2.0,color='k',ls=':')
  axvline(x=sodium_d[0]+linewidth/2.0,color='k',ls=':')
  axvline(x=sodium_d[1]-linewidth/2.0,color='k',ls=':')
  axvline(x=sodium_d[1]+linewidth/2.0,color='k',ls=':')

  xlabel("Wavelength (\AA)")
  xlabel("Wavelength (nm)")
  ylabel('Relative transmission')
  yticks([0.988,0.992,0.996,1.000,1.004,1.008])
  xlim(min_wvl,max_wvl) 
  fig.tight_layout(pad=0.1)
  fig.subplots_adjust(hspace=0)

  tickys = array([5880,5885,5890,5895,5900,5905])

  xticks(tickys)
  ax2.set_xticklabels(tickys/10.0)


  savefig('combined_profile.jpg')
  savefig('combined_profile.pdf')
  os.system('pdftops combined_profile.pdf combined_profile.eps')
  os.system('rm combined_profile.pdf')
  show()

  close()

  fig = plt.figure()

  new_vel = linspace(-500000,500000,116/2)

  vel_specs = []

  spectra = spectra

  vel_resamp_factor = 20


  model = model/baselines

  #imshow(spectra,aspect='auto',cmap='gray',interpolation='none')
  #show()

  #imshow(model,aspect='auto',cmap='gray',interpolation='none')
  #show()

  #imshow((model-spectra),aspect='auto',cmap='gray',interpolation='none')
  #show()

  resamp_factor = 1
  binned_spectra = []
  binned_residuals = []
  binned_model = []
  for i in range(0,len(spectra)):
    bin_wvl, bin_spectra, bin_combined_error = bin_curve(wvl,spectra[i]-1.0,spectra[i],resamp_factor)
    bin_wvl, bin_resid, bin_combined_error = bin_curve(wvl,spectra[i]-model[i],spectra[i],resamp_factor)
    bin_wvl, bin_model, bin_combined_error = bin_curve(wvl,model[i]-1.0,spectra[i],resamp_factor)
    binned_spectra += [super_sample(bin_spectra,resamp_factor)]
    binned_residuals += [super_sample(bin_resid,resamp_factor)]
    binned_model += [super_sample(bin_model,resamp_factor)]

  wvl = wvl[0:len(binned_spectra[0])]

  binned_spectra = array(binned_spectra)
  binned_residuals = array(binned_residuals)
  binned_model = array(binned_model)
  #imshow(binned_spectra,aspect='auto',cmap='gray',interpolation='none',vmin=0.99,vmax=1.01)
  #show()
  #imshow(binned_residuals,aspect='auto',cmap='gray',interpolation='none',vmin=0.99,vmax=1.01)
  #show()

  #spectra = binned_residuals
  spectra = binned_spectra
  #spectra = binned_model

  #for i in range(0,len(spectra)):
    #resamp = rebin_spec(wvl,spectra[i],dummy_wvl)

    #combined = array(list(spectra[i]) + list(spectra[i]))
    #combined_err = array(list(spectra_errors[i]) + list(spectra_errors[i]))
    #combined = combined[vel_index]
    #combined_err = combined_err[vel_index]
    #trial = rebin_spec(combined_vel+10000000,combined,new_vel+10000000)
    #trial_err = rebin_spec(combined_vel+10000000,combined_err,new_vel+10000000)

    #vel_resamp_factor = len(combined_vel[(combined_vel > -500000) & (combined_vel < 500000)])/len(trial)
    #trial = rebin_spec(new_vel+10000000,trial,combined_vel+10000000)
    #trial_err = rebin_spec(new_vel+10000000,trial_err,combined_vel+10000000)

    #bin_vel, bin_combined, bin_combined_error = bin_curve(combined_vel,combined,combined_err,vel_resamp_factor)

    #bin_vel = combined_vel
    #bin_combined = super_sample(bin_combined,vel_resamp_factor)
    #bin_combined_error = super_sample(bin_combined_error,vel_resamp_factor)
    
    #plot(bin_vel,bin_combined)
    #plot(new_vel,trial)
    #show()

    #trial = bin_combined
    #trial_err = bin_combined_error

    #vel_specs += [((trial-1.0)/(trial_err/sqrt(vel_resamp_factor)))]
    #vel_specs += [((trial-1.0)/(trial_err))]
    #vel_specs += [(trial-1.0)]

    #resamp = rebin_spec(wvl,spectra[i],dummy_wvl)
    #resamp_err = rebin_spec(wvl,spectra_errors[i],dummy_wvl)
    #resamp = rebin_spec(wvl,((spectra[i]-1.0)/(spectra_errors[i]/sqrt(resamp_factor))),dummy_wvl)
    #spectra[i] = rebin_spec(dummy_wvl,resamp,wvl)
    #spectra_errors[i] = rebin_spec(dummy_wvl,resamp_err,wvl)
    #spectra[i] = ((spectra[i]-1.0)/(spectra_errors[i]/sqrt(resamp_factor)))

  #new_vel = combined_vel/1000

  #new_vel = bin_vel/1000.0

  #min_vel = -200
  #max_vel = 200

  newspec = []
  #newvelspec = []
  #vel_index = (new_vel > min_vel) & (new_vel < max_vel)

  #print len(vel_index)

  #new_vel = new_vel[vel_index]

  #print len(vel_specs)

  #print len(vel_specs[0])
  #print len(new_vel)
  #print len(vel_index)

  min_wvl = 5887
  max_wvl = 5899
  for i in range(0,len(spectra)):
    index = [(wvl > min_wvl) & (wvl < max_wvl)]
    newspec += [spectra[i][index]]
    #newvelspec += [vel_specs[i][vel_index]]
    #plot(new_vel,newvelspec[-1])
    #show()
  wvl = wvl[index]
  spectra = array(newspec)

  #vel_specs = array(newvelspec)

  #imshow(spectra,aspect='auto',cmap='gray',interpolation='none')

  index = arange(0,len(t))

  #xval = add_2dline(vel_model/1000,new_vel)
  #plot(xval,index,'r--')

  #imshow(vel_specs,aspect='auto',cmap='gray',interpolation='Bicubic')
  #imshow(vel_specs,aspect='auto',cmap='gray',interpolation='none')

  #ax = plt.gca()

  #index = np.arange(99,len(new_vel),200)
  #labels = np.array([int(x) for x in new_vel[index]])
  #xticks(index,labels)

  #t = (time - midtransit)*24
  #index = np.arange(2,len(t),4)
  #labels = np.array([round(x,2) for x in t[index]])
  #yticks(index,labels)

  #ylabel('Time from midtransit (hours)')
  #xlabel("Velocity (km/s)")
  #fig.tight_layout(pad=0.1)

  #cb = colorbar()

  #cb.set_label('Deviation from median')

  #show()

  cmax = amax(binned_spectra)
  cmin= amin(binned_spectra)
  imshow(spectra,aspect='auto',cmap='gray',interpolation='none',vmin=cmin,vmax=cmax)
#  imshow(spectra,aspect='auto',cmap='gray',interpolation='spline36')

  ax = plt.gca()

  index = np.arange(99,len(wvl),200)
  labels = np.array([int(x) for x in wvl[index]])
  xticks(index,labels)

  t = (time - midtransit)*24
  index = np.arange(2,len(t),4)
  labels = np.array([round(x,2) for x in t[index]])
  yticks(index,labels)

  ylabel('Time from midtransit (hours)')
  xlabel("Wavelength (\AA)")
  fig.tight_layout(pad=0.1)

  index = arange(0,len(t))

  line = sin_line(final,time,period,midtransit,wvl,sodium_d[0])
  xval1 = add_2dline(line,wvl)

  line = sin_line(final,time,period,midtransit,wvl,sodium_d[1])
  xval2 = add_2dline(line,wvl)

  cb = colorbar()

  cb.set_label('Deviation from median')

  plot(xval1,index,'r--')
  plot(xval2,index,'r--')

  savefig('2dplot.pdf')
  show()

def add_2dline(xval,xaxis):

  new_xval = xval.copy()
  for i in range(0,len(xval)):
    new_xval[i] = argmin(abs(xval[i]-xaxis))
    print xaxis[new_xval[i]], xval[i]
  return new_xval

def free_func(x,time,period,center,data,errors,wvl,len_pol,plotting=False,linetype='gauss'):

  model = free_model(x,time,period,center,wvl,len_pol,linetype)
  diff = (model-data)/errors

  if plotting == True:
    for i in range(0,len(model)):
      errorbar(wvl,data[i],yerr=errors[i],fmt='bo')
      plot(wvl,model[i],'r-')
      show()

  #imshow(diff)
  #show()

  elements = len(diff)*len(diff[0])

  chi2 = sum(diff**2)/elements

#  print x[0],x[1],chi2

  diff = diff.reshape((elements))
  
  return diff

def voigt_line(x,wvl):
  amp = x[0]
  line_cen = x[1]
  fwhm = 0.4*abs(x[2])
  ratio = abs(x[3])

  Lfwhm = fwhm / (0.5346 + (0.2166 + (ratio**-2))**0.5)
  Gfwhm = fwhm / (0.5346*ratio + (0.2166*ratio**2.0 + 1)**0.5)

  line = voigt(wvl,amp,line_cen,Gfwhm,Lfwhm,normalized=False)

  line = amp*line/max(line)

  return line
  

main()
