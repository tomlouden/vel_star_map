# -*- coding: utf-8 -*-
from pylab import *
import numpy as np
from numpy import pi, cos, sin,sqrt, arcsin, arccos
from spec_manip import *
import scipy.optimize as opt
import pickle
import multiprocessing
import multiprocessing.dummy as multithreading
from functools import partial

def get_star_spec(scale,profile,vsini,u1,u2,spotnumber=0,spotlong=0,spotlat=0,spotsize=0,spot_fill_factor=0,load_star=False,save_star=False,nproc=4,spot_profile=False):

  if load_star == False:
    grid, vel_grid, spots = make_star(scale,vsini,u1,u2,spotnumber,spotlong,spotlat,spotsize,spot_fill_factor)
    star_profile, star_columns, shifted_spectra, shifted_spot_spectra = integrate_star_profile(scale,grid,profile,vel_grid,vsini,spots,nproc=nproc,spot_profile=spot_profile)

    the_star = {'profile':star_profile,'intensity':grid,'velocity':vel_grid,'columns':star_columns,'shifted_spectra':shifted_spectra,'spot_grid':spots,'shifted_spot_spectra':shifted_spot_spectra,}
    if save_star != False:
      pickle.dump(the_star,open(save_star,'wb'))
  else:
    the_star = pickle.load(open(load_star,'rb'))

  return the_star

def vel_prism(data_x, input, profile, times, planet_absorb, spot_data=False, type_data='FLUX',plotting=True,get_star_profile=False,load_star=False,save_star=False,result_wvl=False,nproc=4,location='full',absolute=False,spot_profile=False):

  load_star = False
  # Enter the total stellar and planetary radii  ie (rs + rp) 

  radiusratio = input[0]
  system_scale = input[1]

  semimajor = 1.0/system_scale

  radiustotal = (1.0 + radiusratio)/semimajor

  u1 = input[2]# First LD Coefficent?  

  u2 = input[3]# Second LD Coefficent? 

  inclination = input[4]# Enter the Inclination ie (i)  must be 0 < i > 90
  
  if spot_data == True:
    spotnumber = 1 
  else:
    spotnumber = 0# How many star-spots?  

  # If star spots then create spot propaties with user entered values.

  vsini = input[5]
  planet_K = input[6]
  star_K = input[7]
  midtransit = input[8]
  period = input[9]
  east_offset = input[10]
  west_offset = input[11]
  atm_radius = 1.0 + input[12]

  if spotnumber > 0:
    spotlong = input[13]   # Longituduile postion of spot? ie -90 to 90 note: that 0 is the centre
    spotlat = input[14]  # latitude postion of spot? ie 0 to 90 note: that 90 is the equator
    spotsize = input[15]    # Radius of spot in degrees? note 90 degrees equals the whole stellar disc
    spot_fill_factor = input[16] # Flux of spot? ie 1.00 equals stellar photosphere

  # Convert user variables to system variables

  rs = radiustotal/(radiusratio+1)    # Stellar Radius

  rp = (radiustotal-rs) # gives the planet radii scaled to stellar radius.

  # Now to select a grid size based on system size

  scale = rs/rp

  n = int((100.0 * (rs/rp)))  # Using a planet with a 50 pixel radius 

  # Creates a n*n Grid where n is twice the length of both the planet plus the star radii in pixels

  # Convert inclination to y corrodinates on grid

  inclination_rad = (inclination/180)*(pi) # Converts from degrees to radians

  planet_y_scalar = semimajor * cos(inclination_rad) # finds the planets y corrodinate

  if planet_y_scalar < 0:
    planet_y_scalar = 0 # Cos(90) gives a small negative value, so reset cos(90) = 0

  # set the y-corroditite of the planet

  planet_y = (planet_y_scalar)/(1+(rp/rs)) * ((50.0 * (rs/rp)) +50.0)

  y_offset = (n/2) - planet_y

  # creating the Star on the grid with quadratic limb darkening.

  star_pixel_rad = ((rs/rp) * 50.0)  # Star radius in planetary radii multiplyed by 50 pixels

  x=np.arange(-n/2,n/2)            # Create an vector of length n filled from -n/2 to n/2

  if spotnumber > 0:
    the_star = get_star_spec(scale,profile,vsini,u1,u2,spotnumber=1,spotlong=spotlong,spotlat=spotlat,spotsize=spotsize,spot_fill_factor=spot_fill_factor,load_star=load_star,save_star=save_star,nproc=nproc,spot_profile=spot_profile)
  else:
    the_star = get_star_spec(scale,profile,vsini,u1,u2,load_star=load_star,save_star=save_star,nproc=nproc,spot_profile=spot_profile)


  vel_grid = the_star['velocity']
  grid = the_star['intensity']
  star_profile = the_star['profile']
  star_columns = the_star['columns']
  shifted_spectra = the_star['shifted_spectra']

  wvl = profile['wvl']

  vel_model = planet_K*sin(2.0*pi*(midtransit-times)/period)
  star_vel_model = -1*star_K*sin(2.0*pi*(midtransit-times)/period)


  # Intergrations of the planet to create LC model.

#  let's multiprocess this instead, set it up as a new function and map k to pool tada?...

  p = multiprocessing.Pool(nproc)

  data_points = len(data_x)

  # make curry

  curry = []
  for i in range(0,data_points):
    curry += [[n,data_points,data_x,semimajor,star_pixel_rad,rs,rp,atm_radius,planet_y,grid,the_star,profile,x,times,vel_model,star_vel_model,east_offset,west_offset,planet_absorb,star_columns,star_profile,shifted_spectra,result_wvl,plotting,1,i,location,absolute,spot_profile]]

  if nproc == 1:
    output = []
    for i in range(0,data_points):
      results = uncurry_model_transit(curry[i])
      output += [results]
  else:
    output = p.map(uncurry_model_transit,curry)
    
    # Release memory
    p.close()

  output = array(output)

  return output

def make_star(scale,vsini,u1,u2,spotnumber=0,spotlong=0,spotlat=0,spotsize=0,spot_fill_factor=0):

  star_pixel_rad = (scale * 50.0)  # Star radius in planetary radii multiplyed by 50 pixels
  n = int((100.0 * scale))  # Using a planet with a 50 pixel radius 

  grid = np.zeros((n,n))
  vel_grid = np.zeros((n,n))

  x=np.arange(-n/2,n/2)            # Create an vector of length n filled from -n/2 to n/2
  for i in range(0,n-1):          # Begin small filling loop
    r = sqrt(x**2+(i-(n/2.0))**2)     # creates an x by i array filled with r 
    larray = np.arange(0,len(x))
    star = larray[r <= star_pixel_rad]  # Creates star array. if r > star_pixel_rad returns a single value of -1. Else creates a vector with all the x values at slice i which fall inside the star.  

    if(len(star) > 0):
      m = cos(arcsin(r[star]/star_pixel_rad))
      grid[star,i]= (1-u1*(1-m)-u2*(1-m)**2)  # Fills the grid slice i at the x values generated by star with the star and the appropirate LD value.
#      grid[star,i]= 1.0  # Fills the grid slice i at the x values generated by star with the star and the appropirate LD value.

  vels = vsini*cos(pi*1.0*arange(0,n-1)/(n-1))

  for i in range(0,n-1):          # Begin small filling loop
    r = sqrt(x**2+(i-(n/2.0))**2)     # creates an x by i array filled with r 
    larray = np.arange(0,len(x))
    star = larray[r <= star_pixel_rad]  # Creates star array. if r > star_pixel_rad returns a single value of -1. Else creates a vector with all the x values at slice i which fall inside the star.  

    if(len(star) > 0):
      m = cos(arcsin(r[star]/star_pixel_rad))
      vel_grid[i,star]= vels[i]  # Fills the grid slice i at the x values generated by star with the star and the appropirate LD value.

  #Star Spots

  spot = np.zeros((n,n))
  if spotnumber > 0: 

    for sp in range(0,len(spotlong)):
        plane_centre = np.zeros(3)
        # need to find x and y corordinates of spots
        spotlong_rad = (pi * (spotlong[sp]) / 180)
        spotlat_rad = (pi * (spotlat[sp]) / 180)
        spotsize_rad = (pi * (spotsize[sp]) / 180)
        spx = star_pixel_rad * sin(spotlong_rad) * sin(spotlat_rad)
        spy = star_pixel_rad * cos(spotlat_rad)
        spz = star_pixel_rad * cos(spotlong_rad) * sin(spotlat_rad)
        sps = star_pixel_rad * sin(spotsize_rad)
        spot_centre = np.array([spx,spy,spz])
        plane_centre = spot_centre * cos(spotsize_rad)
    
  ##############################Creating the Starspot#############################################################

        for i in range(int(spx-1.1*sps), int(spx + 1.1*sps)):
          for j in range(int(spy-1.1*sps), int(spy+1.1*sps)):
            r = sqrt((i-spx)**2 + (j-spy)**2)
            if r < star_pixel_rad:
              xxx = i
              yyy = j
              zzz = sqrt((star_pixel_rad)**2 - xxx**2 - yyy**2)
              point_centre = np.array([xxx,yyy,zzz])
              diff_vector = np.zeros(3)
              diff_vector = point_centre - plane_centre
              dot_product = sum(plane_centre * diff_vector)
              diff_mag = sqrt(diff_vector[0]**2 + diff_vector[1]**2 + diff_vector[2]**2)
              plane_mag = sqrt(plane_centre[0]**2 + plane_centre[1]**2 + plane_centre[2]**2)
              theta = arccos(dot_product / (diff_mag * plane_mag))
              if ((theta <= pi/2.0) and (i+n/2 <= n-1) and (j+n/2 <= n-1)):
                spot[i+n/2,j+n/2] = spot_fill_factor[sp]

#        grid = grid * spot

  return grid, vel_grid, spot


def uncurry_model_transit(c):
  return model_transit(c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12],c[13],c[14],c[15],c[16],c[17],c[18],c[19],c[20],c[21],c[22],c[23],c[24],c[25],c[26],c[27],c[28])

def model_transit(n,data_points,data_x,semimajor,star_pixel_rad,rs,rp,atm_radius,planet_y,grid,the_star,profile,oldx,times,vel_model,star_vel_model,east_offset,west_offset,planet_absorb,star_columns,star_profile,shifted_spectra,result_wvl,plotting,nproc,k,location='full',absolute=False,spot_profile=False):

      planet = np.zeros((n,n)) + 1
      atmosphere = np.zeros((n,n)) + 1

      new_data_x = data_x.copy()

      # Fully optermized method using strips to scan over 100 * 100 pixel box and only for data points In transit

      velocities = []
      actual_vel = []
      time_so_far = []

      xx = ( (sin(new_data_x[k] * 2 * pi)*(semimajor * star_pixel_rad))) +((50.0 * (rs/rp)))
      y = np.arange(100*atm_radius) + planet_y - 50*atm_radius

      for i in range(int(xx - 51),int(xx + 51)):
          r = sqrt((y - planet_y)**2 + (i-xx)**2)
          index = arange(0,len(r))
          planet_disc = index[r <= 50.0]
          if (len(planet_disc) > 0) and (i > 0) and (i < n-1):
              index = ((planet_disc + ((50.0 * (rs/rp))) + planet_y) - 50.0*atm_radius)
              index = [int(x) for x in index]
              index = array(index)
              index = index[index < n]
              planet[i,index] = 0

      for i in range(int(xx - 51*atm_radius), int(xx + 51*atm_radius)):
          r = sqrt((y - planet_y)**2 + (i-xx)**2)
          index = arange(0,len(r))
          planet_disc = index[(r > 50.0) & (r <= 50.0*atm_radius)]
          if (len(planet_disc) > 0) and (i > 0) and (i < n-1):
              index = ((planet_disc + ((50.0 * (rs/rp))) + planet_y) - 50.0*atm_radius)
              index = [int(x) for x in index]
              index = array(index)
              index = index[index < n]
              if i < xx:
                atmosphere[i,index] = 0.5
              else:
                atmosphere[i,index] = -0.5

      

      wvl = profile['wvl']
      flux = grid * planet    
      new_flux = sum(sum(flux))

      east_planet_absorb = redshift(wvl,planet_absorb,(vel_model[k] + east_offset)*1e3)
      west_planet_absorb = redshift(wvl,planet_absorb,(vel_model[k] + west_offset)*1e3)

      # MIGHT WANT TO DO SOMETHING ABOUT STELLAR ORB VELOCITY AT SOME POINT, BUT NOT YET
      #shifted_spectra = redshift(wvl,shifted_spectra,star_vel_model[k]*1e3)


      offset = (east_offset + west_offset)/2
      ang_vel = (east_offset - offset)

      if plotting == True:
        import matplotlib.pyplot as plt
        contrast = sum(spot_profile)/sum(profile['spectrum'])
        picture = flux*(1.0 - the_star['spot_grid']) + contrast*flux*(the_star['spot_grid']) 
        plt.imshow(picture.T,cmap='gist_heat')
        plt.savefig('figure.pdf')
        plt.close()

      r_planet_profile = integrate_rotation_profile(profile,flux,the_star,shifted_spectra,oldx,n,star_pixel_rad,planet_absorb,xx,atm_radius,offset,ang_vel,vel_model[k],atmosphere=True,atmosphere_grid=atmosphere,star_columns=star_columns,nproc=nproc,location=location)

      if absolute == False:
        r_differential_transit = (r_planet_profile/median(r_planet_profile) - (star_profile/median(star_profile)))/(star_profile/median(star_profile))
        differential_transit = r_differential_transit
        if any(result_wvl != False):
          results = rebin_spec(wvl,differential_transit + 1,result_wvl)
        else:
          results = differential_transit + 1

      else:
        r_differential_transit = r_planet_profile/median(star_profile) - star_profile/median(star_profile)
        differential_transit = r_differential_transit
        if any(result_wvl != False):
          results = rebin_spec(wvl,differential_transit + 1,result_wvl)
        else:
          results = differential_transit + 1

#      differential_transit = (planet_profile/median(planet_profile) - (star_profile/median(star_profile)))/(star_profile/median(star_profile))

#      plot(differential_transit)
#      plot(r_differential_transit)
#      show()
#          quit()

      results = r_planet_profile/star_profile

      return results

def illum_track(data_x, input, profile, times, planet_absorb, spot_data=False, type_data='FLUX',plotting=True,get_star_profile=False,load_star=False,save_star=False,result_wvl=False,nproc=4,location='full',absolute=False):

  load_star = False

  radiusratio = input[0]
  system_scale = input[1]

  semimajor = 1.0/system_scale

  radiustotal = (1.0 + radiusratio)/semimajor

  u1 = input[2]# First LD Coefficent?  

  u2 = input[3]# Second LD Coefficent? 

  inclination = input[4]# Enter the Inclination ie (i)  must be 0 < i > 90
  
  if spot_data == True:
    spotnumber = 1 
  else:
    spotnumber = 0# How many star-spots?  

  # If star spots then create spot propaties with user entered values.

  vsini = input[5]
  planet_K = input[6]
  star_K = input[7]
  midtransit = input[8]
  period = input[9]
  east_offset = input[10]
  west_offset = input[11]
  atm_radius = 1.0 + input[12]

  if spotnumber > 0:
    spotlong = input[13]   # Longituduile postion of spot? ie -90 to 90 note: that 0 is the centre
    spotlat = input[14]    # latitude postion of spot? ie 0 to 90 note: that 90 is the equator
    spotsize = input[15]    # Radius of spot in degrees? note 90 degrees equals the whole stellar disc
    spot_fill_factor = input[16] # Flux of spot? ie 1.00 equals stellar photosphere

  # Convert user variables to system variables

  rs = radiustotal/(radiusratio+1)    # Stellar Radius

  rp = (radiustotal-rs) # gives the planet radii scaled to stellar radius.

  # Now to select a grid size based on system size

  scale = rs/rp

  n = int((100.0 * (rs/rp)))  # Using a planet with a 50 pixel radius 

  # Creates a n*n Grid where n is twice the length of both the planet plus the star radii in pixels

  # Convert inclination to y corrodinates on grid

  inclination_rad = (inclination/180)*(pi) # Converts from degrees to radians

  planet_y_scalar = semimajor * cos(inclination_rad) # finds the planets y corrodinate

  if planet_y_scalar < 0:
    planet_y_scalar = 0 # Cos(90) gives a small negative value, so reset cos(90) = 0

  # set the y-corroditite of the planet

  planet_y = (planet_y_scalar)/(1+(rp/rs)) * ((50.0 * (rs/rp)) +50.0)

  y_offset = (n/2) - planet_y

  # creating the Star on the grid with quadratic limb darkening.

  star_pixel_rad = ((rs/rp) * 50.0)  # Star radius in planetary radii multiplyed by 50 pixels

  x=np.arange(-n/2,n/2)            # Create an vector of length n filled from -n/2 to n/2

  if spotnumber > 0:
    the_star = get_star_spec(scale,profile,vsini,u1,u2,spotnumber=0,spotlong=0,spotlat=0,spotsize=0,load_star=load_star,save_star=save_star,nproc=nproc)
  else:
    the_star = get_star_spec(scale,profile,vsini,u1,u2,load_star=load_star,save_star=save_star,nproc=nproc)


  vel_grid = the_star['velocity']
  grid = the_star['intensity']
  star_profile = the_star['profile']
  star_columns = the_star['columns']
  shifted_spectra = the_star['shifted_spectra']

  wvl = profile['wvl']

  star_rv = 2.1

  vel_model = star_rv + 1*planet_K*sin(2.0*pi*(midtransit-times)/period)
  star_vel_model = -1*star_K*sin(2.0*pi*(midtransit-times)/period)


  # Intergrations of the planet to create LC model.

#  let's multiprocess this instead, set it up as a new function and map k to pool tada?...

  p = multiprocessing.Pool(nproc)

  data_points = len(data_x)

  # make curry

  new_data_x = data_x.copy()

  xx = ( (sin(new_data_x * 2 * pi)*(semimajor * star_pixel_rad))) +((50.0 * (rs/rp)))
  old_xx = xx.copy()
  thickness = 50*atm_radius
  if location == 'left':
    xx = xx - thickness
  if location == 'right':
    xx = xx + thickness


  int_xx = np.array([int(x) for x in xx])
  int_y = int(planet_y)

  d1 = 5895.924
  d2 = 5889.950

  c = 299792.458

  shifted_d1 = d1*(1.0 - (vel_model/c))
  shifted_d2 = d2*(1.0 - (vel_model/c))

  background_d1 = []
  background_d2 = []
  star_vel = []
  out_shifted_spectra = []

  for i in range(0,len(vel_model)):
    if (xx[i] > max(old_xx)) or (xx[i] < min(old_xx)):
      background_d1 += [0]
      background_d2 += [0]
      star_vel += [0]
      out_shifted_spectra += [shifted_spectra[0].copy()*0.0]
    else:
      ld = grid[xx[i],planet_y]
      med_shifted_spectra = ld*shifted_spectra[int_xx[i]]/np.median(shifted_spectra[0])
      background_d1 += [med_shifted_spectra[argmin(abs(wvl-shifted_d1[i]))]]
      background_d2 += [med_shifted_spectra[argmin(abs(wvl-shifted_d2[i]))]]
      star_vel += [vel_grid[int_xx[i],planet_y]]
      out_shifted_spectra += [shifted_spectra[int_xx[i]]]


  background_d2 = np.array(background_d2)
  background_d1 = np.array(background_d1)
  star_vel = np.array(star_vel)

  return wvl, background_d1, background_d2, star_vel, vel_model-star_rv, shifted_d1,out_shifted_spectra

def integrate_star_profile(scale,grid,profile,vel_grid,vsini,spots,nproc=4,spot_profile=False):
  
  spectrum = profile['spectrum'].copy()
  wvl = profile['wvl'].copy()

  star_pixel_rad = (scale * 50.0)  # Star radius in planetary radii multiplyed by 50 pixels
  n = int((100.0 * scale))  # Using a planet with a 50 pixel radius 
  x=np.arange(-n/2,n/2)            # Create an vector of length n filled from -n/2 to n/2

  curry = []
  for i in range(0,n-1):
    curry += [[i,x,n,star_pixel_rad,grid,spectrum,wvl,vel_grid,vsini,spots,spot_profile]]

  if nproc > 1:
    p = multiprocessing.Pool(nproc)
    output = p.map(uncurry_profile_multi,curry)
    p.close()
  else:
    output = []
    for i in range(0,n-1):
      output += [uncurry_profile_multi(curry[i])]

  output = array(output)
  combined_profile = sum(output[:,0],axis=0)
  star_columns = output[:,0]
  shifted_spectra = output[:,1]
  shifted_spot_spectra = output[:,2]

#  return combined_profile/median(combined_profile), star_columns, shifted_spectra
  return combined_profile, star_columns, shifted_spectra, shifted_spot_spectra

def integrate_rotation_profile(profile,grid,the_star,shifted_spectra,x,n,star_pixel_rad,planet_absorb,xx,atm_radius,offset,ang_vel,k_vel,atmosphere=False,atmosphere_grid=False,star_columns=[],nproc=4,plotting=False,location='full'):

  combined_profile = profile['spectrum'].copy()*0.0
  spectrum = profile['spectrum'].copy()
  wvl = profile['wvl'].copy()

  in_atmosphere = []
  not_in_atmosphere = []
  Longitude = []

  vel_grid = the_star['velocity']
  spots = the_star['spot_grid']
  shifted_spot_spectra = the_star['shifted_spot_spectra']
  shifted_star_spectra = the_star['shifted_spectra']
    
  for i in range(0,n-1):          # Begin small filling loop
    r = sqrt(x**2+(i-(n/2.0))**2)     # creates an x by i array filled with r 
    larray = np.arange(0,len(x))
    star = larray[r <= star_pixel_rad]  # Creates star array. if r > star_pixel_rad returns a single value of -1. Else creates a vector with all the x values at slice i which fall inside the star.  

    if(len(star) > 0):
      shifted = profile['spectrum'].copy()*0.0
      corrected = profile['spectrum'].copy()*0.0

      if any(abs(atmosphere_grid[i,star]) == 0.5):
        vel_model = k_vel + offset + ang_vel*sin(0.5*pi*(xx-i)/(50*atm_radius))
        if len(wvl) > 100:
          shifted_planet_absorb = redshift(wvl,planet_absorb,vel_model*1e3)
        else:
          shifted_planet_absorb = planet_absorb
        Longitude += [(xx-i)/(50*atm_radius)]
        for j in range(0,len(star)):

          if spots[i,star[j]] == 0:
            shifted_spectra = shifted_star_spectra[i]
          else:
            shifted_spectra = spots[i,star[j]]*shifted_spot_spectra[i] +  (1 - spots[i,star[j]])*shifted_star_spectra[i]

          if abs(atmosphere_grid[i,star[j]]) == 0.5:
            shifted += shifted_spectra*grid[i,star[j]]
          else:
            corrected += shifted_spectra*grid[i,star[j]]

        absorbed = shifted*shifted_planet_absorb
        in_atmosphere += [absorbed]
        not_in_atmosphere += [shifted]

      elif any(grid[i,star] == 0):
        for j in range(0,len(star)):

          if spots[i,star[j]] == 0:
            shifted_spectra = shifted_star_spectra[i]
          else:
            shifted_spectra = spots[i,star[j]]*shifted_spot_spectra[i] +  (1 - spots[i,star[j]])*shifted_star_spectra[i]

          corrected += shifted_spectra*grid[i,star[j]]

      else:
        corrected = star_columns[i]

      combined_profile += corrected


  in_atmosphere = np.array(in_atmosphere)
  not_in_atmosphere = np.array(not_in_atmosphere)
  Longitude = np.array(Longitude)

  summed_atmosphere = combined_profile.copy()
  summed_atmosphere = 0.0*summed_atmosphere

  for ii in range(0,len(Longitude)):
    if location == 'leftest':
      if Longitude[ii] == min(Longitude):
        summed_atmosphere += in_atmosphere[ii]*(len(Longitude))
      else:
        summed_atmosphere += not_in_atmosphere[ii]

  for ii in range(0,len(Longitude)):
    if location == 'left':
      if Longitude[ii] < 0:
        summed_atmosphere += in_atmosphere[ii]
      else:
        summed_atmosphere += not_in_atmosphere[ii]


    if location == 'rightest':
      if Longitude[ii] == max(Longitude):
        summed_atmosphere += in_atmosphere[ii]*(len(Longitude))
      else:
        summed_atmosphere += not_in_atmosphere[ii]

    if location == 'right':
      if Longitude[ii] > 0:
        summed_atmosphere += in_atmosphere[ii]
      else:
        summed_atmosphere += not_in_atmosphere[ii]

    if location == 'full':
      summed_atmosphere += in_atmosphere[ii]

    if location == 'none':
      summed_atmosphere += not_in_atmosphere[ii]



  combined_profile += summed_atmosphere

  return combined_profile


def uncurry_profile_multi(c):
  return star_profile_multi(c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10])

def star_profile_multi(i,x,n,star_pixel_rad,grid,spectrum,wvl,vel_grid,vsini,spots,spot_profile):

    r = sqrt(x**2+(i-(n/2.0))**2)     # creates an x by i array filled with r 
    larray = np.arange(0,len(x))
    star = larray[r <= star_pixel_rad]  # Creates star array. if r > star_pixel_rad returns a single value of -1. Else creates a vector with all the x values at slice i which fall inside the star.  
    corrected = spectrum.copy()*0.0
    if vsini != 0:
      shift_spec = redshift(wvl,spectrum,median(vel_grid[i,star])*1e3)
      shift_spot_spec = redshift(wvl,spot_profile,median(vel_grid[i,star])*1e3)
    else:
      shift_spec = spectrum
      shift_spot_spec = spot_profile
    for j in range(0,len(star)):
      if spots[i,star[j]] == 0.0:
        corrected += shift_spec*grid[i,star[j]]
      else:
        corrected += spots[i,star[j]]*shift_spot_spec*grid[i,star[j]] +  (1 - spots[i,star[j]])*shift_spec*grid[i,star[j]]
    return corrected, shift_spec, shift_spot_spec