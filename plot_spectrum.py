# -*- coding: utf-8 -*-
import pyfits as pf
import os
from pylab import *
from scipy.optimize import leastsq
from bcar import boxcar
import matplotlib.pyplot as plt
from numpy import *

def main():

  # selected wavelength ranges

  fold_factor = 10

#  name = 'corot-2b'
#  name = 'corot-3b'
#  name = 'corot-18b'
  name = 'HD189733b'
#  name = 'WASP-2b'
#  name = 'WASP-4b'
#  name = 'WASP-5b'
#  name = 'WASP-6b'
#  name = 'WASP-7b'
#  name = 'WASP-8b'
#  name = 'WASP-15b'
#  name = 'WASP-16b'
#  name = 'WASP-17b'
#  name = 'WASP-18b'
#  name = 'WASP-19b'
#  name = 'WASP-22b'
#  name = 'WASP-23b'
#  name = 'WASP-24b'
#  name = 'WASP-25b'
#  name = 'WASP-26b'
#  name = 'WASP-30b'
#  name = 'WASP-31b'
#  name = 'WASP-32b'
#  name = 'WASP-38b'
#  name = 'WASP-71b'
#  name = 'WASP-80b'

  data_vals = []

  for line in open('param_files/'+name):
    data_vals += [line.split('=')[1].strip('\n').strip(' ')]

  Target = data_vals[0].split('.or.')
  period = float(data_vals[1])
  tmid = float(data_vals[2])
  duration = float(data_vals[3])*24.0
  date_query = data_vals[4]

  data_dir = './data/'+Target[0]+'/'

  filelist = os.listdir(data_dir)
  filelist = [data_dir+x for x in filelist]

  s1d_list =[x for x in filelist if 's1d' in x and 'HARPS' in x]
  s1d_list =[x for x in filelist if 's1d' in x]
  s1d_list_A =[x for x in s1d_list if '_A' in x and date_query in x]

  tmid -= 2400000.5

  wvl_range = []

  # sodium
  wvl_range += [[5887,5899]]
  # comparison
  wvl_range += [[5875,5887]]
  wvl_range += [[5899,5911]]

# white
#  wvl_range += [[4500,6500]]
  


  wvl_range = array(wvl_range)

  s1d_list_B =[x for x in s1d_list if '_B' in x and '2007-08-' in x]

#  s1d_list_A =[x for x in s1d_list if '_A' in x]

#  s1d_list_A =[x for x in s1d_list_A if '2006-07' in x or '2006-08' in x or '2006-09' in x or '2007-08' in x]
#  s1d_list_A =[x for x in s1d_list_A if '2006-07' in x or '2006-08' in x]

#  ccf_check(s1d_list_A,Target,tmid,period)
#  quit()

  create_difference_image(s1d_list_A,Target,tmid,period,duration,fold_factor)

#  quit()

  mjd, total_flux, normed_flux = extract_data(s1d_list_A,wvl_range,Target)

  transits = ephem_calc([min(mjd)-1.0,max(mjd)+1.0],tmid,period)

  hour = 1.0/24.0

  for transit in transits:
    depth = [0,100]
    plot([transit,transit],depth,'k-')
    plot([transit+hour*duration/2.0,transit+hour*duration/2.0],depth,'k:')
    plot([transit-hour*duration/2.0,transit-hour*duration/2.0],depth,'k:')

  for i in range(0,len(wvl_range)):
#    plot(mjd,total_flux[:,i]/median(total_flux[:,i]),'o')	
    plot(mjd,normed_flux[:,i],'o')

  ylabel('FLUX (Relatvie CPS)')
  xlabel('MJD')

  show()

  sodium_transit = ((normed_flux[:,0]) - (normed_flux[:,2] + normed_flux[:,1]))/(normed_flux[:,2] + normed_flux[:,1])

  for transit in transits:
    depth = [-0.639,-0.645]
    plot([transit,transit],depth,'k-')
    plot([transit+hour*duration/2.0,transit+hour*duration/2.0],depth,'k:')
    plot([transit-hour*duration/2.0,transit-hour*duration/2.0],depth,'k:')

  ylabel('Sodium - White light / white (Relatvie CPS)')
  plot(mjd,sodium_transit,'o')
  show()

#  mjd_2, total_flux_2, normed_flux_2 = extract_data(s1d_list_B,wvl_range,ref)
#  for i in range(0,len(wvl_range)):
#    plot(mjd,total_flux[:,i]/median(total_flux[:,i]),'o')	
#    plot(mjd_2,median(np.append(normed_flux_2[:10,i],normed_flux_2[-8:,i])) - normed_flux_2[:,i],'x')

def extract_data(file_list,wvl_range,Target):

  mjd = []
  total_flux = []
  normed_flux = []

  for filen in file_list:
    with pf.open(filen) as HDU:
      header = HDU[0].header 
      data = HDU[0].data
# work out conversion factor
      ref = header['CRVAL1']
      inc = header['CDELT1']

#      wvl = ref + (np.arange(0,len(data))*inc)
#      plot(wvl,data)
#      show()
#      quit()

      exptime = header['EXPTIME']

      print filen, header['OBJECT']
      range_fluxes = [] 
      normed_fluxes = [] 

      if any([x in header['OBJECT'] for x in Target]):    
        mjd += [header['MJD-OBS']]
        for i in wvl_range:
          pix_range = [int(round(x)) for x in (i - ref)/inc]
          range_fluxes += [sum(data[pix_range[0]:pix_range[1]])/exptime]
          normed_fluxes += [sum(data[pix_range[0]:pix_range[1]])/((pix_range[1] - pix_range[0])*exptime)]
        total_flux += [range_fluxes]
        normed_flux += [normed_fluxes]
      else:
	print 'this is not the star you are looking for', header['OBJECT']

  return array(mjd), array(total_flux), array(normed_flux)

def ccf_check(file_list,Target,tmid,period):

  rvs = []
  phase = []
  for filen in file_list:
    print filen
    try:
      with pf.open(filen.replace('s1d','ccf_G2')) as HDU:
	if any([x in HDU[0].header['OBJECT'] for x in Target]):
	  rvs += [fit_ccf(HDU)]
	  phase += [find_phase(HDU,tmid,period)]
    except:
      with pf.open(filen.replace('s1d','ccf_K5')) as HDU:
	if any([x in HDU[0].header['OBJECT'] for x in Target]):
	  rvs += [fit_ccf(HDU)]
	  phase += [find_phase(HDU,tmid,period)]

  print 'using '+str(len(rvs))+'rv measurements'

  rvs = array(rvs)
  phase = array(phase)

  zero_point = rvs[argmin(abs(phase))]

  plot(phase,(rvs - zero_point)*1000,'rx')
  ylabel('Radial Velocity (m/s)')
  xlabel('Phase')
  show()

def find_phase(HDU,tmid,period):

  mjd = HDU[0].header['MJD-OBS']
  
  nearest_transit = ephem_calc([mjd-period/2.0,mjd+period/2.0],tmid,period)

  return (mjd - nearest_transit[0])/period

def fit_ccf(image):
#  ccf_min = argmin(image[0].data[-1])
  ccf = image[0].data[-1]
# prior = [half way along ccf, height of ccf, 1/10th length of ccf, median of ccf]
  prior = [len(ccf)/2.0,max(ccf) - min(ccf), len(ccf)/10.0, median(ccf), 0.01, 0.01]
  x, cov = leastsq(gauss_absorb, prior, args=(ccf))

#  gaussian_line = gaussian(x,ccf)
#  plot(ccf - gaussian_line,'rx')
#  plot(gaussian_line - gaussian_line)
#  show()

  ccf_min = x[0]
  ref = image[0].header['CRVAL1']
  ind = image[0].header['CDELT1']
  rv = ccf_min*ind + ref
  return rv

def gauss_absorb(x,data):
  gaussian_line = gaussian(x,data)
  return gaussian_line - data

def gaussian(x,data):
  dummy = np.arange(0,len(data))
  return x[1]*exp(-((dummy - x[0])**2)/(2*x[2]**2)) + x[3] + x[4]*dummy + x[5]*(dummy)**2.0

def create_difference_image(file_list,Target,tmid,period,duration,fold_factor):

  with pf.open(file_list[0]) as HDU:
    header = HDU[0].header 
    data = HDU[0].data
# work out conversion factor
    ref = header['CRVAL1']
    inc = header['CDELT1']

    wvl = ref + (np.arange(0,len(data))*inc)

  inc_temp = inc
  ref_temp = ref

  in_transit = array([0.0]*len(wvl))
  pre_transit = array([0.0]*len(wvl))
  post_transit = array([0.0]*len(wvl))

  n = [0,0]

  for filen in file_list:
    with pf.open(filen) as HDU:
      header = HDU[0].header 
      data = array(HDU[0].data)
      #check if is transit
      date = header['MJD-OBS']
      search_time = [date-period/2.0,date+period/2.0]
      transit = ephem_calc(search_time,tmid,period)
      if((date < transit + (duration/24.0)/2.0) and (date > transit - (duration/24.0)/2.0)) and any([x in header['OBJECT'] for x in Target]):
#      if((date > transit + (duration/24.0)/2.0) and Target in header['OBJECT']):
	shift = int(round((header['CRVAL1'] - ref_temp)/inc_temp))
	n[0] += 1.0
	print 'in transit', filen
        if shift > 0:
	  data = np.append([0]*abs(shift),data)
        if shift < 0:
	  data = data[abs(shift):]
        data = array(np.append(data,array([0]*abs(len(data) - len(in_transit)))))
#	data = data/header['EXPTIME']
	in_transit += data[0:len(in_transit)]
	plot(wvl,data[0:len(in_transit)]/sum(data[0:len(in_transit)]),'g')
	ylim([-0.000001,0.000006])
	xlim([5888,5892])
	savefig('plot_bin/'+filen.split('/')[-1]+'.png', bbox_inches=0)
	clf()


      else:
	if any([x in header['OBJECT'] for x in Target]):
	  shift = int(round((header['CRVAL1'] - ref_temp)/inc_temp))
	  print 'not in transit', filen
	  n[1] += 1.0
	  if shift > 0:
	    data = np.append([0]*abs(shift),data)
	  if shift < 0:
	    data = data[abs(shift):]
	  data = array(np.append(data,array([0]*abs(len(data) - len(pre_transit)))))
#	data = data/header['EXPTIME']
	  if((date < transit - (duration/24.0)/2.0) and any([x in header['OBJECT'] for x in Target])):
	    pre_transit += data[0:len(pre_transit)]
	    plot(wvl,data[0:len(pre_transit)]/sum(data[0:len(pre_transit)]),'b')
	    ylim([-0.000001,0.000006])
	    xlim([5888,5892])
	    savefig('plot_bin/'+filen.split('/')[-1]+'.png', bbox_inches=0)
	    clf()
	  else:
	    post_transit += data[0:len(post_transit)]
	    plot(wvl,data[0:len(post_transit)]/sum(data[0:len(post_transit)]),'b')
	    ylim([-0.000001,0.000006])
	    xlim([5888,5892])
	    savefig('plot_bin/'+filen.split('/')[-1]+'.png', bbox_inches=0)
	    clf()


  pre_transit = array([float(x) for x in pre_transit])
  post_transit = array([float(x) for x in post_transit])
  in_transit = array([float(x) for x in in_transit])/sum(in_transit)

  out_transit = (pre_transit + post_transit) / sum((pre_transit + post_transit))
  pre_transit = pre_transit / sum(pre_transit)
  post_transit = post_transit / sum(post_transit)

  print n

  out_transit = out_transit
  in_transit = in_transit

  error_out_transit = sqrt(out_transit)
  error_in_transit = sqrt(in_transit)

  difference = (in_transit - out_transit)/out_transit

  difference_null = (pre_transit - post_transit)/post_transit

  difference_error = (((sqrt(error_out_transit**2.0 + error_in_transit**2.0))/(out_transit + in_transit))**2.0 + (error_out_transit/out_transit)**2.0)**0.5

  print difference_error

#  plot(difference_error)

  xlimit = [5860,5930]
  ylimit = [-0.2,0.2]

  sodium_d = [5889.950,5895.924]

  plot(wvl,out_transit)
  plot(wvl,in_transit)

  ylabel("Signal (Total 'relative' counts)")
  xlabel("Wavelength (Angstroms)")

  xlim(xlimit)
  show()

  box_wvl, box_difference = boxcar(wvl,difference,fold_factor)

  box_wvl_null, box_difference_null = boxcar(wvl,difference_null,fold_factor)

  plot(wvl,difference)
  ylabel("IN - OUT / OUT ('Relative' counts per second)")
  xlabel('Wavelength (Angstroms)')
  xlim(xlimit)
  ylim(ylimit)
  show()

  
  fig, (ax0,ax1) = plt.subplots(nrows=2, sharex=True)

  ax0.plot(box_wvl,box_difference)
  ax0.set_ylabel("IN - OUT / OUT")
  ax0.set_xlabel('Wavelength (Angstroms)')
  ax0.set_ylim(array(ylimit)/10.0)
  ax0.plot([sodium_d[0],sodium_d[0]],[-1,1],'k:')
  ax0.plot([sodium_d[1],sodium_d[1]],[-1,1],'k:')


  ax1.plot(box_wvl_null,box_difference_null)
  ax1.set_ylabel("PRE - POST / POST")
  ax1.set_xlabel('Wavelength (Angstroms)')
  ax1.set_xlim(xlimit)
  ax1.set_ylim(array(ylimit)/10.0)
  ax1.plot([sodium_d[0],sodium_d[0]],[-1,1],'k:')
  ax1.plot([sodium_d[1],sodium_d[1]],[-1,1],'k:')

  show()

  bin_factor = 33

  bin_wvl, bin_difference = bin_curve(wvl,difference,bin_factor)
  bin_wvl, bin_difference_error = bin_curve(wvl,difference_error,bin_factor)

  bin_difference_error = bin_difference_error/sqrt(bin_factor)

  white = array(bin_difference[4500/bin_factor:6850/bin_factor])

  white = white[white != 'nan']
 
  print 'The mean of the difference is: ', mean(white)

  errorbar(bin_wvl,bin_difference,yerr=bin_difference*(bin_difference_error),fmt='bD')

  ylabel("IN - OUT / OUT ('Relative' counts per second)")
  xlabel("Wavelength (Angstroms)")

  xlim(xlimit)
  ylim(ylimit)
  show()

def bin_curve(axis,curve,curve_error,factor,avgtype='mean'):
  nbins = int(len(curve)/factor)
  change = int(len(curve) - nbins*factor)
  oldcurve = curve.copy()-1.0

  if change != 0:
    axis=axis[:-change]
    curve=curve[:-change]
    curve_error=curve_error[:-change]
  axis =array([mean(axis[x:x+factor]) for x in xrange(1, len(axis), factor)])

  if avgtype == 'mean':
    curve =array([mean(curve[x:x+factor]) for x in xrange(1, len(curve), factor)])
  if avgtype == 'median':
    curve =array([median(curve[x:x+factor]) for x in xrange(1, len(curve), factor)])

  #two different ways of estimating errors
  curve_error =array([mean(curve_error[x:x+factor]) for x in xrange(1, len(curve_error), factor)])
  curve_error = curve_error/sqrt(1.0*factor)

  #curve_error =array([std(oldcurve[x:x+factor]) for x in xrange(1, len(curve_error), factor)])
  #curve_error = curve_error/sqrt(1.0*factor)

  return axis, curve, curve_error

def ephem_calc(search_time,tmid,period):
  search_grid = tmid + period*np.arange(-1000,1000)
  return search_grid[(search_grid > search_time[0]) & (search_grid < search_time[1])]


#main()
