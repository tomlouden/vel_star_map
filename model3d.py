# -*- coding: utf-8 -*-
from numpy import loadtxt
from numpy import array
from numpy import sin
from numpy import pi
#from pyspeckit.spectrum.models.inherited_voigtfitter import voigt
import numpy as np
from vel_prism import vel_prism, illum_track
import scipy

def fit_model_3d(x,time,period,planet_K,star_K,midtransit,spectra,spectra_errors,wvl,line_centers,plotting=False,nproc=4,best_fit=False):


  if plotting == True:
    save_star = 'star_bin.p'
  else:
    save_star = False
  
  model = model_3d(x,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,nproc=nproc,best_fit=best_fit,save_star=save_star)
  diff = (model-spectra)/spectra_errors

  vel_model = planet_K*sin(2.0*pi*(midtransit-time)/period)

  d1 = 5895.924
  d2 = 5889.950
  velocity = 10e3
  c = 3e8

  resamp_factor = 5
  if plotting == True:
    x_no_atm = x.copy()
    best_fit['atm_radius'] = 0.0
    no_atm_model = model_3d(x_no_atm,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,nproc=nproc,best_fit=best_fit,load_star=save_star)

    for i in range(0,len(model)):
      bin_wvl, bin_spectra, bin_combined_error = bin_curve(wvl,spectra[i],spectra_errors[i],resamp_factor)
      errorbar(bin_wvl,bin_spectra,yerr=bin_combined_error,fmt='k.')
      plot(wvl,model[i],'r-')
      plot(wvl,model[i]/no_atm_model[i],'g-')
      xlim(5887,5899)
      ylim(0.97,1.03)
      savefig('west_atm_plots/'+str(i)+'.png')
      clf()

    #errorbar(wvl,mean((spectra/no_atm_model),axis=0),std((spectra/no_atm_model),axis=0)/sqrt(len(spectra)),fmt='k.')
    #plot(wvl,mean((model/no_atm_model),axis=0),'r-')
    #show()

    alligned = align_spectra(wvl,spectra.copy(),vel_model)
    alligned_model = align_spectra(wvl,model.copy(),vel_model)
    corrected_model = align_spectra(wvl,model.copy()/no_atm_model,vel_model)

    bin_wvl, bin_spectra, bin_combined_error = bin_curve(wvl,mean(spectra,axis=0),std(spectra,axis=0)/sqrt(len(spectra)),resamp_factor)
    errorbar(bin_wvl,bin_spectra,yerr=bin_combined_error,fmt='k.')
    plot(wvl[10:-10],mean(alligned_model,axis=0)[10:-10],'r-')
    plot(wvl[10:-10],mean(corrected_model,axis=0)[10:-10],'g-')
    show()

  #imshow(diff)
  #show()

  elements = len(diff)*len(diff[0])

  chi2 = sum(diff**2)/elements

#  print x[0],x[1],chi2

  print(x)
  print(chi2)

  diff = diff.reshape((elements))
  
  return diff

def model_3d(x,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,star_profile,planet_absorb,nproc=4,best_fit=False,save_star=False,load_star=False,star_vsini=0,location='full',absolute=False,result_wvl=False,spot_profile=False):


    p0 = 0.1572
    system_scale = 1.0/8.92
    gamma_0 = x[0]
    gamma_1 = x[1]
    inc = 85.68
    east_offset = x[2]*0
    west_offset = x[3]*0
    atm_radius = x[4]
    fwhm = x[5]
    ratio = x[6]
    atm_strength = x[7]



#    if len(x) == 8:
#      p0 = 0.1572
#      system_scale = 1.0/8.92
#      gamma_0 = x[0]
#      gamma_1 = x[1]
#      inc = 85.68
#      east_offset = x[2]*0
#      west_offset = x[3]*0
#      atm_radius = x[4]
#      fwhm = x[5]
#      ratio = x[6]
#      atm_strength = x[7]

#    else:
#      p0 = 0.1572
#      system_scale = 1.0/8.92
#      gamma_0 = best_fit['ld1']
#      gamma_1 = best_fit['ld2']
#      inc = 85.68
      #star_vsini = 3.1
#      east_offset = x[0]
#      west_offset = x[1]

#      atm_radius = best_fit['atm_radius']
#      fwhm = x[2]
#      ratio = x[3]
#      atm_strength = x[4]

#      gamma_0 = x[0]
#      gamma_1 = x[1]
#      east_offset = x[2]*0
#      west_offset = x[3]*0
#      atm_radius = best_fit['atm_radius']
#      fwhm = x[4]
#      ratio = x[5]
#      atm_strength = x[6]


    spotlong, spotlat, spotsize, spotflux = x[8],x[9],x[10],x[11]

    spotlong = []
    spotlat = []
    spotsize = []
    spotflux = []
    for i in range(0,len(x[8:])/4):
      spotlong +=[x[i*4 + 8]]
      spotlat +=[x[i*4 + 9]]
      spotsize +=[x[i*4 + 10]]
      spotflux +=[x[i*4 + 11]]

    input = [p0,system_scale,gamma_0,gamma_1,inc,star_vsini,planet_K,star_K,midtransit,period,east_offset,west_offset,atm_radius,spotlong, spotlat, spotsize, spotflux]
    data_x = (time-midtransit)/period
#    master_dat = loadtxt('sodium_spectrum.dat')
#    model_wvl = array(master_dat[:,0])
#    master_flux = array(master_dat[:,1])
#    profile = {'wvl':model_wvl,'spectrum':master_flux}

    #plot(model_wvl,planet_absorb)
    #show()
    #quit()
    spot_data = True

    results = vel_prism(data_x, input, star_profile, time,planet_absorb,spot_data=spot_data, type_data='FLUX',plotting=plotting,result_wvl=result_wvl,nproc=nproc,save_star=save_star,load_star=load_star,location=location,absolute=absolute,spot_profile=spot_profile)
    return results

def illumination_model(x,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,star_profile,nproc=4,best_fit=False,save_star=False,load_star=False,star_vsini=3.1,location='full',absolute=False):

    if len(x) == 8:
      p0 = 0.1572
      system_scale = 1.0/8.92
      gamma_0 = x[0]
      gamma_1 = x[1]
      inc = 85.68
      #star_vsini = 3.1
      east_offset = x[2]
      west_offset = x[3]
      atm_radius = x[4]
      fwhm = x[5]
      ratio = x[6]
      atm_strength = x[7]

    else:
      p0 = 0.1572
      system_scale = 1.0/8.92
#      gamma_0 = best_fit['ld1']
#      gamma_1 = best_fit['ld2']
      inc = 85.68
      #star_vsini = 3.1
      east_offset = x[0]
      west_offset = x[1]

      atm_radius = best_fit['atm_radius']
      fwhm = x[2]
      ratio = x[3]
      atm_strength = x[4]

      gamma_0 = x[0]
      gamma_1 = x[1]
      east_offset = x[2]
      west_offset = x[3]
      atm_radius = best_fit['atm_radius']
      fwhm = x[4]
      ratio = x[5]
      atm_strength = x[6]


    input = [p0,system_scale,gamma_0,gamma_1,inc,star_vsini,planet_K,star_K,midtransit,period,east_offset,west_offset,atm_radius]
    data_x = (time-midtransit)/period
#    master_dat = loadtxt('sodium_spectrum.dat')
#    model_wvl = array(master_dat[:,0])
#    master_flux = array(master_dat[:,1])
#    profile = {'wvl':model_wvl,'spectrum':master_flux}

    planet_absorb = array([1.0]*len(model_wvl))
    amp = 1.0
    fwhm = 0.4*abs(fwhm)
    ratio = abs(ratio)
    strengths = [1.0,0.5]

    s = 0
    for center in line_centers:
      planet_absorb = planet_absorb*(1.0 - atm_strength*strengths[s]*voigt_line([amp,center,fwhm,ratio],model_wvl))
      s += 1

    for i in range(0,len(planet_absorb)):
      if planet_absorb[i] < 0:
        planet_absorb[i] = 0

    #plot(model_wvl,planet_absorb)
    #show()
    #quit()

    results = illum_track(data_x, input, profile, time,planet_absorb,spot_data=False, type_data='FLUX',plotting=plotting,result_wvl=wvl,nproc=nproc,save_star=save_star,load_star=load_star,location=location,absolute=absolute)
    return results


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

def voigt(xarr,amp,xcen,sigma,gamma,normalized=False):
    """
    Normalized Voigt profile

    z = (x+i*gam)/(sig*sqrt(2))
    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))

    The area of V in this definition is 1.
    If normalized=False, then you can divide the integral of V by
    sigma*sqrt(2*pi) to get the area.

    Original implementation converted from 
    http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
    (had an incorrect normalization and strange treatment of the input
    parameters)

    Modified implementation taken from wikipedia, using the definition.
    http://en.wikipedia.org/wiki/Voigt_profile

    Parameters
    ----------
    xarr : np.ndarray
        The X values over which to compute the Voigt profile
    amp : float
        Amplitude of the voigt profile
        if normalized = True, amp is the AREA
    xcen : float
        The X-offset of the profile
    sigma : float
        The width / sigma parameter of the Gaussian distribution
    gamma : float
        The width / shape parameter of the Lorentzian distribution
    normalized : bool
        Determines whether "amp" refers to the area or the peak
        of the voigt profile
    """

    z = ((xarr-xcen) + 1j*gamma) / (sigma * np.sqrt(2))
    V = amp * np.real(scipy.special.wofz(z)) 
    if normalized:
        return V / (sigma*np.sqrt(2*np.pi))
    else:
        return V
        #tmp = 1.0/scipy.special.wofz(numpy.zeros((len(xarr))) \
        #      +1j*numpy.sqrt(numpy.log(2.0))*Lfwhm).real
        #tmp = tmp*amp* \
        #      scipy.special.wofz(2*numpy.sqrt(numpy.log(2.0))*(xarr-xcen)/Gfwhm+1j* \
        #      numpy.sqrt(numpy.log(2.0))*Lfwhm).real
        #return tmp

