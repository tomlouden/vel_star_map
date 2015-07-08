# -*- coding: utf-8 -*-
from model3d import *
import os
import pickle
from pylab import *
import numpy as np
from plot_spectrum import bin_curve
import scipy.signal as signal
import matplotlib

def main():


  x1 = [0.51374991, 0.18869174, 7.1395660212000003, 0.34985208915100002, 0.23553298503299999, 0.39732301091299999, 9.2615599123899998e+44, 5.1084213145500001]
  x1 = [0.51374991, 0.18869174, 8.5049162065200008, -0.15442891679199999,  0.23553298503299999, 0.063145416345300007, -3.50757872893e+50, 111.002838225]
  label1 = 'ld real' #standard errors

  #x1 = [0.51374991, 0.18869174, 7.5354499164600002, -0.095576620693199998, 0.22283249803800001, 0.54247514116999995, 1.20815021896e+45, 3.5895451398099998]
  #label1 = 'ld real' #small errors

  x2 = [2.16157235,-1.72168638,6.94563836,-0.02163093,0.23282932,0.53659169,9.87500318e+44,4.09058055]
  x2 = [2.3793539848950003, -2.3551000484399998, 7.6581408896499994, 1.0522893867450001, 0.23553298503299999, 0.28236854826800001, -3.6695715340100001e+50, 9.5370933805650004]
  x2 = [2.3793539848950003, -2.3551000484399998, 7.6581408896499994, 1.0522893867450001, 0, 0.28236854826800001, -3.6695715340100001e+50, 9.5370933805650004]
  x2 = [0.51374991, 0.18869174, 4.17524364486, 4.17524364486,  0.23553298503299999, 0.063145416345300007, -3.50757872893e+50, 111.002838225]
  label2 = 'ld free'
  #x2 = [1.76975025,-0.76728519,7.66635193,-0.09572382,0.21485677,0.58513335,1.21818912e+45,3.66011491]
  #label2 = 'ld fixed'

  outputf = 'model.avi'
  outputf2 = 'periodigram.avi'

  #x1 = [2.16157235,-1.72168638,6.94563836,-0.02163093,0.23282932,0.53659169,9.87500318e+44,4.09058055]
  #label1 = 'rotation'
  #x2 = [2.16157235,-1.72168638,3.462003715,3.462003715,0.23282932,0.53659169,9.87500318e+44,4.09058055]
  #label2 = 'constant offset'
  #outputf = 'compare_rotation.avi'

  dir = 'everything/'

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
    index = [(data['wvl'] > min_wvl) & (data['wvl'] < max_wvl)]

    wvl += [data['wvl'][index]]
    spectra += [data['spec'][index]]
    spectra_errors += [data['spec_error'][index]]

  wvl = wvl[0]
  time = array(time)
  spectra = array(spectra)
  spectra_errors = array(spectra_errors)


  planet_K = 154
  sodium_d = [5889.950,5895.924]
  midtransit = 54341.056842
  period =  2.21857312
  planet_K = 154
  star_K = 0.20196
  line_centers = sodium_d
  plotting = False
  nproc = 4
  best_fit=False
  save_star = False

  i1 = 9
  i2 = 29


  spectra = spectra[argsort(time)]
  spectra = spectra[np.arange(i1,i2)]

  spectra_errors = spectra_errors[argsort(time)]
  spectra_errors = spectra_errors[np.arange(i1,i2)]

  time = time[argsort(time)]
  time = time[np.arange(i1,i2)]

  new_errors = array([std(spectra, axis=0)]*20)

  print spectra_errors
  print new_errors

  spectra_errors = new_errors

  load = True
  save = True
  save_m1 = 'model1.p'
  save_m2 = 'model2.p'

  if load == False:
    #model1 = model_3d(x1,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,nproc=nproc,best_fit=best_fit,save_star=save_star,star_vsini=3.1)
    #model2 = model_3d(x2,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,nproc=nproc,best_fit=best_fit,save_star=save_star,star_vsini=3.1)

    model1 = model_3d(x1,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,nproc=nproc,best_fit=best_fit,save_star=save_star,star_vsini=0)
    model2 = model_3d(x2,time,period,planet_K,star_K,midtransit,wvl,line_centers,plotting,nproc=nproc,best_fit=best_fit,save_star=save_star,star_vsini=0)

    if save == True:
      pickle.dump(model1, open(save_m1, 'wb'))
      pickle.dump(model2, open(save_m2, 'wb'))
      
  else:
    model1 = pickle.load(open(save_m1))
    model2 = pickle.load(open(save_m2))


  os.system('rm model_plot/*.png')

  #spectra_errors = array(spectra_errors)*sqrt(0.5)

  #for i in range(0,len(model1)):
    #plot(wvl,model1[i],'r-')
    #plot(wvl,model2[i],'g-')
    #show()

  flist = ''
  flist2 = ''
  delta_chi2 = []


  bad_chi2 = sum(((spectra-1)/spectra_errors)**2,axis=0)
  good_chi2 = sum(((spectra-model1)/spectra_errors)**2,axis=0)


  ffour, axone = plt.subplots(2, 1,sharex=True)

  good_scaled_chi = sum(((spectra-model1)/spectra_errors),axis=0)

  axone[0].errorbar(wvl,good_scaled_chi/sqrt(20),yerr=1,fmt='k,',alpha=1.0,capsize=0,elinewidth=0.5,marker=None)

  #axone[1].errorbar(wvl,good_chi2-bad_chi2,yerr=1,fmt='k,',alpha=1.0,capsize=0,elinewidth=0.5,marker=None)
  axone[1].plot(wvl,good_chi2-bad_chi2,'k-')

  axone[1].set_xlim(5888.0,5898.0)

  axone[0].set_ylim(-13,13)
  axone[1].set_ylim(-23,13)

  xlabel('Wavelength ($\AA$)')

  axone[0].set_ylabel('Normalised residuals')
  axone[1].set_ylabel('$\Delta {\chi}^2$')

  #axone[0].yaxis.set_ticks([])

  gcf().subplots_adjust(bottom=0.15)

  ffour.subplots_adjust(hspace=0,wspace=0)

  savefig('deltachi_comparison.png')

  savefig('deltachi_comparison.pdf')

  close()

  bad_absolute = sum(abs((spectra-1)/spectra_errors),axis=0)
  good_absolute = sum(abs((spectra-model1)/spectra_errors),axis=0)

  bad_chi = sum(((spectra-1)),axis=0)
  good_chi = sum(((spectra-model1)),axis=0)

  bad_scaled_chi = sum(((spectra-1)/spectra_errors),axis=0)
  good_scaled_chi = sum(((spectra-model1)/spectra_errors),axis=0)

  #ffour, axfour = plt.subplots(3, 3,sharex=True)

  #axfour = axfour.reshape(-1)


  #axfour[0].plot(wvl,bad_absolute)
  #axfour[3].plot(wvl,good_absolute)
  #axfour[6].plot(wvl,good_absolute-bad_absolute)

  #axfour[1].plot(wvl,bad_chi/sqrt(20))
  #axfour[4].plot(wvl,good_chi/sqrt(20))
  #axfour[7].plot(wvl,(bad_chi-good_chi)/sqrt(20))

  #axfour[2].plot(wvl,bad_scaled_chi/sqrt(20))
  #axfour[5].plot(wvl,good_scaled_chi/sqrt(20))
  #axfour[8].plot(wvl,(bad_scaled_chi-good_scaled_chi)/sqrt(20))

  #axfour[4].set_xlim(5880.0,5905.0)
  
  #show()

  for plotn in range(0,2):

    fbig, axbig = plt.subplots(5, 4)
    fbig.set_size_inches(8.27,11.69)

    fbig.subplots_adjust(hspace=0.1,wspace=0.1)

    axbig = axbig.reshape(-1)

    resamp_factor = 6


    residuals = []
    null_residuals = []
    

    for i in range(0,len(model1)/2):

      si = i +(10*plotn)

      binned_wvl, binned_spectra, binned_spectra_errors = bin_curve(wvl,spectra[si],spectra_errors[si],resamp_factor)

      chiresid = (spectra[si]-model1[si])/spectra_errors[si]

      residuals += [(spectra[si]-model1[si])]
      null_residuals += [(spectra[si]-(model1[si]/model1[si]))]

      axbig[i*2].set_xlim(5888.0,5892.0)
      axbig[i*2].errorbar(wvl,spectra[si],yerr=spectra_errors[si],fmt='k,',alpha=1.0,capsize=0,elinewidth=0.5,marker=None)
      #axbig[i*2].errorbar(binned_wvl,binned_spectra,yerr=binned_spectra_errors,fmt='k.',capsize=0)
      axbig[i*2].plot(wvl,model1[si],'r')
      axbig[i*2].set_ylim(0.940,1.060)
      axbig[i*2].set_ylim(0.880,1.060)
      axbig2 = axbig[i*2].twinx()
      #axbig2.plot(wvl,chiresid,'k,')
      axbig2.errorbar(wvl,chiresid,yerr=chiresid/chiresid,fmt='k,',alpha=1.0,capsize=0,elinewidth=0.5,marker=None)
      axbig2.plot(wvl,-1 + (model1[si]/model1[si]),'r')
      axbig2.set_ylim(-5,20)
      axbig[i*2].set_xlim(5888.0,5892.0)
      axbig2.set_xlim(5888.0,5892.0)

      axbig[(i*2)+1].set_xlim(5894.0,5898.0)
      axbig[(i*2)+1].errorbar(wvl,spectra[si],yerr=spectra_errors[si],fmt='k,',alpha=1.0,capsize=0,elinewidth=0.5,marker=None)
      #axbig[(i*2)+1].errorbar(binned_wvl,binned_spectra,yerr=binned_spectra_errors,fmt='k.',capsize=0)
      axbig[(i*2)+1].plot(wvl,model1[si],'r')
      axbig[(i*2)+1].set_ylim(0.880,1.060)
      axbig3 = axbig[(i*2)+1].twinx()
      #axbig3.plot(wvl,chiresid,'k,')
      axbig3.errorbar(wvl,chiresid,yerr=chiresid/chiresid,fmt='k,',alpha=1.0,capsize=0,elinewidth=0.5,marker=None)
      axbig3.plot(wvl,-1 + (model1[si]/model1[si]),'r')
      axbig3.set_ylim(-5,20)
      axbig[(i*2)+1].set_xlim(5894.0,5898.0)
      axbig3.set_xlim(5894.0,5898.0)

      axbig[i*2].spines['right'].set_visible(False)
      axbig[i*2].yaxis.tick_left()

      axbig2.tick_params(labelright='off')
      axbig2.set_yticks([])
      axbig[(i*2)+1].set_yticks([])

      if i % 2 != 0:
	axbig[i*2].tick_params(labelright='off')
	axbig[i*2].tick_params(labelleft='off')
	if i == 5:
	  axbig3.set_ylabel('Normalised residuals',rotation=270,labelpad=15)
      if i % 2 == 0:
	if i == 4:
	  axbig[(i*2)].set_ylabel('Differential spectra')
	axbig2.tick_params(labelleft='off')
	axbig3.tick_params(labelright='off')
	axbig3.tick_params(labelleft='off')

      range1 = [5888,5889,5890,5891]
      range2 = [5894,5895,5896,5897]


      axbig2.xaxis.set_ticks(range1)
      axbig3.xaxis.set_ticks(range2)
      axbig[i*2].xaxis.set_ticks(range1)
      axbig[(i*2)+1].xaxis.set_ticks(range2)

      if (10 - i) > 2:
	axbig2.xaxis.set_ticklabels([])
	axbig3.xaxis.set_ticklabels([])
	axbig[i*2].xaxis.set_ticklabels([])
	axbig[(i*2)+1].xaxis.set_ticklabels([])
      else:
	axbig2.xaxis.set_ticklabels(range1)
	axbig3.xaxis.set_ticklabels(range2)
	axbig[i*2].xaxis.set_ticklabels(range1)
	axbig[(i*2)+1].xaxis.set_ticklabels(range2)
	axbig[i*2].set_xlabel('          Wavelength ($\AA$)')
	axbig[i*2].xaxis.set_label_coords(1,-0.11)

      axbig2.spines['right'].set_visible(False)
      axbig[(i*2)+1].spines['left'].set_visible(False)
      axbig[(i*2)+1].tick_params(labelleft='off')
      axbig3.spines['left'].set_visible(False)

      for item in (axbig2.get_xticklabels() + axbig2.get_yticklabels() + axbig3.get_xticklabels() + axbig3.get_yticklabels() + axbig[(i*2)+1].get_xticklabels() + axbig[(i*2)+1].get_yticklabels() + axbig[i*2].get_xticklabels() + axbig[i*2].get_yticklabels()):
	item.set_fontsize(12)


      d = .03

      kwargs = dict(transform=axbig[(i*2)].transAxes, color='k', clip_on=False)
      axbig[(i*2)].plot((1-d,1+d),(-d,+d), **kwargs)
      axbig[(i*2)].plot((1-d,1+d),(1-d,1+d), **kwargs)

      kwargs = dict(transform=axbig[(i*2)+1].transAxes, color='k', clip_on=False)
      axbig[(i*2)+1].plot((-d,+d),(-d,+d), **kwargs)
      axbig[(i*2)+1].plot((-d,+d),(1-d,1+d), **kwargs)

    savefig('bigfig_'+str(plotn)+'.pdf')
    savefig('bigfig_'+str(plotn)+'.png')

    close()


  # stop here for now
  return

  sum_null = mean(null_residuals, axis=0)/(sqrt(3)*spectra_errors[si]/sqrt(si))
  sum_residuals = mean(residuals, axis=0)/(sqrt(3)*spectra_errors[si]/sqrt(si))

  sum_null = sum((null_residuals/spectra_errors[si])**2, axis=0)
  sum_residuals = sum((residuals/spectra_errors[si])**2, axis=0)

  sum_null = sum(null_residuals/spectra_errors[si], axis=0)
  sum_residuals = sum(residuals/spectra_errors[si], axis=0)

  #sum_null = mean(null_residuals, axis=0)
  #sum_residuals = mean(residuals, axis=0)

  #errorbar(wvl,sum_null,yerr=spectra_errors[si]/spectra_errors[si],fmt='k,',alpha=1.0,capsize=0,elinewidth=0.5,marker=None)

  #stack = mean(((spectra[:10])),axis=0)
  #stack_model = mean(((model1[:10])),axis=0)
  #binned_wvl, binned_spectra, binned_spectra_errors = bin_curve(wvl,stack,spectra_errors[i]/sqrt(10/3),5)
  #errorbar(binned_wvl,binned_spectra,yerr=binned_spectra_errors,fmt='ko')
  #plot(wvl,stack_model,'r-')

  #stack = mean(((spectra[10:])),axis=0)
  #stack_model = mean(((model1[10:])),axis=0)
  #binned_wvl, binned_spectra, binned_spectra_errors = bin_curve(wvl,stack,spectra_errors[i]/sqrt(10/3),5)
  #errorbar(binned_wvl,binned_spectra-0.03,yerr=binned_spectra_errors,fmt='ko')
  #plot(wvl,stack_model-0.03,'r-')

  show()
  
  for i in range(0,len(model1)):

    ni = str(i)
    if len(ni) < 2:
      ni = '0'+str(i)

    chi2diffmap = ((model1[i] - spectra[i])/spectra_errors[i])**2 - ((model2[i] - spectra[i])/spectra_errors[i])**2
    #plot(wvl,np.cumsum(chi2diffmap),label='difference')
    #close()
  
    chi2_1 = sum(((model1[i] - spectra[i])/spectra_errors[i])**2)
    chi2_2 = sum(((model2[i] - spectra[i])/spectra_errors[i])**2)
    delta_chi2 += [chi2_1 - chi2_2]
    print label1, label2
    print chi2_1, chi2_2
    print chi2_1/len(model1[i]), chi2_2/len(model1[i])
    print delta_chi2[-1]


    chiresid = (spectra[i]-model1[i])/spectra_errors[i]

    plot_freq =False
    if plot_freq == True:
      freqs = linspace(0.01,20,1000)
      x = arange(0,len(model1[i]))
      pgram = signal.lombscargle(x.astype('float64'),chiresid.astype('float64'),freqs.astype('float64'))
      f2 , ax5 = plt.subplots(1, 1)
      ax5.plot(freqs,pgram)
      fname = 'model_plot/'+ni+'_periodigram.png'
      savefig(fname)
      close()
      flist2 += fname+','

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    binned_wvl, binned_spectra, binned_spectra_errors = bin_curve(wvl,spectra[i],spectra_errors[i],resamp_factor)


    for ax in [ax1, ax2]:
      ax.errorbar(wvl,spectra[i],yerr=spectra_errors[i],fmt='k,',alpha=0.2,capsize=0)
      ax.errorbar(binned_wvl,binned_spectra,yerr=binned_spectra_errors,fmt='k.',capsize=0)
      ax.plot(wvl,model1[i],'r',label=label1)
      #ax.plot(wvl,model2[i],'g--',label=label2)
      ax.set_ylim(0.940,1.070)

    for ax in [ax3, ax4]:  
      ax.errorbar(wvl,chiresid,yerr=spectra_errors[i]/spectra_errors[i],fmt='k,',alpha=0.2,capsize=0)
      binned_wvl, binned_spectra, binned_spectra_errors = bin_curve(wvl,(spectra[i]-model1[i])/spectra_errors[i],spectra_errors[i]/spectra_errors[i],resamp_factor)
      ax.errorbar(binned_wvl,binned_spectra,yerr=binned_spectra_errors,fmt='k.',capsize=0)
      ax.plot(wvl,-1 + model1[i]/model1[i],'r',label=label1)
      ax.set_ylim(-5,5)

    ax1.set_xlabel('Wavelength ($\AA$)')
    ax1.set_xlim(5880,5905)
    ax1.set_xlim(5887.5,5897.5)

    ax1.set_xlim(5888.0,5892.0)
    ax3.set_xlim(5888.0,5892.0)

    ax2.set_xlim(5894,5898)
    ax4.set_xlim(5894,5898)

    f.suptitle(str(round(delta_chi2[-1],2))+' delta chi2 '+' rchi1 '+str(round(chi2_1/len(model1[i]),2))+' rchi2 '+str(round(chi2_2/len(model1[i]),2)))
    #legend()

    f.suptitle('rchi2 '+str(round(chi2_1/len(model1[i]),2)))
    gcf().subplots_adjust(bottom=0.15)
    fname = 'model_plot/'+ni+'.png'
    savefig(fname)
    flist += fname+','

    close()
  flist = flist.rstrip(',')

  print sum(delta_chi2), 'total delta chi2'

  print flist
  fps = 1

  command = 'mencoder mf://'+flist+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
  print command
  os.system(command)
  os.system('firefox '+ outputf)

  command = 'mencoder mf://'+flist2+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf2
  print command
  os.system(command)
  #os.system('firefox '+ outputf2)

main()