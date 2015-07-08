# -*- coding: utf-8 -*-
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
import sys
from scipy.stats import ttest_ind
from scipy.stats import ks_2samp
from numpy import *
import numpy as np
import os
from make_transit_diagram import add_diagram



rcParams.update({'font.size': 12})

def main():

  linetype = 'lorentz'

  infile = sys.argv[-1]

  if linetype == 'voigt':
    absolute_list = [False,False,False,True,True]
    plot_lims = [(-50,250),(-15,5),(0,0.015),(0,1.0),(-10,10)]
    for i in range(0,5):
      make_plots(i,plot_lims,infile,absolute=absolute_list[i])
  else:
    absolute_list = [False,False,False,True,True]
    plot_lims = [(-50,250),(-10,10),(0,0.015),(0,1.0),(-10,10)]
    for i in range(1,2):
      make_plots(i,plot_lims,infile,absolute=absolute_list[i])

def make_plots(par_num,plot_lims,infile,absolute=False):

  fig = plt.figure()
  my_locator = MaxNLocator(6)
  plt.minorticks_on()

  out_names = ['K_P','V_P','Depth','FWHM','Voigt_param']
  par_names = ['$K_P$ (km s$^{-1}$)','$V_P$ (km s$^{-1}$)','Fractional depth','FWHM (\AA)','voigt param (dimensionless)']

  par_label = par_names[par_num]
  out_name = out_names[par_num]

  i = 0

  data_dict = {}

  burns = 0 

  for line in open(infile):
    if i ==0:
      data_names = [str(x) for x in line.strip(' \n').split(' ')]
      for d in data_names:
	data_dict[d] = []
    else:
      if i > burns:
	d = [float(x) for x in line.strip(' \n').split(' ')]
	for x in range(0,len(d)):
	  data_dict[data_names[x]] += [d[x]]
    i += 1

  #data = loadtxt(openf)

  print len(data_dict['eastvel']), 'samples'

  try:
    atm = array(data_dict['atm'])
  except:
    atm = 0.235380980973

  rotation = 2.6
  factor = 2.0/pi

  print factor

  westvel = ((array(data_dict['eastvel'])- 2.3)*-1)
  eastvel = ((array(data_dict['westvel']) - 2.3)*-1)
  avgvel = (eastvel + westvel)/2

  data_dict['eastvel'] = (1.0 + (atm/2.0))*array(data_dict['eastvel'])/(1.0 + atm)
  data_dict['westvel'] = (1.0 + (atm/2.0))*array(data_dict['westvel'])/(1.0 + atm)

  westvel = ((array(data_dict['eastvel'])- 2.3)*-1)
  eastvel = ((array(data_dict['westvel']) - 2.3)*-1)

  diff = abs(eastvel-westvel)

  percent = np.percentile(eastvel, [0.15,2.5,16, 50, 84,97.5,99.85],axis=0)
  percent2 = np.percentile(westvel+2.91475245023, [0.15,2.5,16, 50, 84,97.5,99.85],axis=0)
  percent3 = np.percentile(avgvel, [0.15,2.5,16, 50, 84,97.5,99.85],axis=0)

  percent4 = np.percentile(diff, [0.15,2.5,16, 50, 84,97.5,99.85],axis=0)

  try:
    atm = array(data_dict['atm'])
  except:
    atm = 0.235380980973

  period = 2.2185733

  RJ=69911.0

  prad = 1.138*RJ

  aprad = prad*(1.0 + (atm/2.0))

  apcircum = aprad*2*pi

  rotvel = apcircum/(period*24*3600)

  print rotvel

  windvel = diff - (rotvel*2)

  percent5 = np.percentile(windvel, [0.15,2.5,16, 50, 84,97.5,99.85],axis=0)

  print 'eastvel'
  print percent
  print 'westvel'
  print percent2
  print 'avgvel'
  print percent3
  print 'diffvel'
  print percent4
  print 'windvel'
  print percent5


  print len(diff[diff > 0])/(1.0*len(diff)), len(diff[diff < 0])/(1.0*len(diff))

  #hist(diff)
  #show()

  print 'from ',len(eastvel),'trials'

  bins = 20

  axp = subplot(2,1,1)

  axp.axis('off')

  u1 = 0.51374991000000003
  u2 = 0.18869174

  #u1 = 2.4
  #u2 = -2.46

  add_diagram(axp,u1,u2)

  ax1 = subplot(2,3,4)
  ax2 = subplot(2,3,5)
  ax3 = subplot(2,3,6)

  ax1.hist(eastvel,histtype='stepfilled',color=(1,0,0),alpha=0.1,bins=bins,edgecolor=(1,0,0))
  ax1.hist(eastvel,histtype='step',color=(1,0,0),alpha=1,bins=bins,edgecolor=(1,0,0))
  #ax1.hist(param_1,histtype='step',color=(1,0,0),alpha=1,bins=bins,hatch='////',edgecolor=(1,0,0))

  ax2.hist(avgvel,histtype='stepfilled',color=(0.5,0,0.5),alpha=0.1,bins=bins,edgecolor=(0.5,0,0.5))
  ax2.hist(avgvel,histtype='step',color=(0.5,0,0.5),alpha=1,bins=bins,edgecolor=(0.5,0,0.5))
  #ax2.hist(param_1_3,histtype='step',color=(0.5,0,0.5),alpha=1,bins=bins,hatch='|||',edgecolor=(1,0,1))

  ax3.hist(westvel,histtype='stepfilled',color=(0,0,1),alpha=0.1,bins=bins,edgecolor=(0,0,1))
  ax3.hist(westvel,histtype='step',color=(0,0,1),alpha=1,bins=bins,edgecolor=(0,0,1))
  #ax3.hist(param_1_2,histtype='step',color=(0,0,1),alpha=1,bins=bins,hatch='\\\\\\\\',edgecolor=(0,0,1))

  xlabel_list = [-20,-10,0,10,20]

  xlabel_list = [-10,-5,0,5]

  ylabel_list = [1000,2000,3000]

  ax1.set_ylabel('Number')
  ax1.set_xticks(xlabel_list)
  ax1.set_xticklabels(xlabel_list)

  #line = ax1.get_lines()[0]
  #y0 = min(line.get_ydata())
  #y1 = max(line.get_ydata())

  #print y1, y0

  #quit()

  ax2.set_xlabel(par_label)
  #ax2.set_ylabel('Number')
  ax2.yaxis.set_ticklabels([])
  ax2.set_xticks(xlabel_list)
  ax2.set_xticklabels(xlabel_list)

  #ax3.set_ylabel('Number')
  ax3.yaxis.set_ticklabels([])
  ax3.set_xticks(xlabel_list)
  ax3.set_xticklabels(xlabel_list)

  ax1.set_xlim(plot_lims[par_num])
  ax2.set_xlim(plot_lims[par_num])
  ax3.set_xlim(plot_lims[par_num])

  #red_patch = mpatches.Patch(color=(1,0,0),alpha=0.1)
  red_step = mpatches.Patch(ec=(1,0,0),fill=False,hatch='////')

  #blue_patch = mpatches.Patch(color=(0,0,1),alpha=0.1)
  blue_step = mpatches.Patch(ec=(0,0,1),fill=False,hatch='\\\\\\\\')

  green_step = mpatches.Patch(ec=(0,1,0),fill=False,hatch='|||')

  #lg = plt.legend([blue_step,green_step,red_step],['Egress','Full transit','Ingress'])
  #lg.draw_frame(False)

  fig.tight_layout(pad=0.1)

  lw = 2
  best_fit = -116.566826293
  #for p in percent:
    #axvline(x=p,color='k',ls='--',linewidth=lw)
    #axvline(x=p,color='k',ls='--',linewidth=lw)
  #axvline(x=best_fit,color='r',ls='-',linewidth=lw)

  ax1.axvline(x=0,color='k',ls='--')
  ax2.axvline(x=0,color='k',ls='--')
  ax3.axvline(x=0,color='k',ls='--')

  plt.gcf().subplots_adjust(bottom=0.15)

  out_name = out_name + '_' + infile.rstrip('.txt')
  print out_name
  savefig(out_name+'.png')
  savefig(out_name+'.pdf')
  os.system('pdftops '+out_name+'.pdf '+out_name+'.eps')

main()
