# -*- coding: utf-8 -*-
from pylab import *
import triangle
import sys

nwalkers = 16

#openf = 'test_mcmc.txt'
#openf = 'mcmc_no_bootstrap.txt'

openf = sys.argv[-2]
burns = float(sys.argv[-1])

#openf = 'correct_ld.txt'


data = []
i = 0

data_dict = {}

for line in open(openf):
  if i ==0:
    data_names = [str(x) for x in line.strip(' \n').split(' ')]
    for d in data_names:
      data_dict[d] = []
  else:
    if i > burns:
      d = [float(x) for x in line.strip(' \n').split(' ')]
      data += [d]
      for x in range(0,len(d)):
        data_dict[data_names[x]] += [d[x]]
  i += 1

data = array(data)

print shape(data)

#data = loadtxt(openf)

print len(data_dict['eastvel']), 'samples'
print len(data)/nwalkers, 'iterations'

try:
  atm = array(data_dict['atm'])
except:
  atm = 0.235380980973

perc = percentile(data,[2.275,15.865,50,84.135,97.725],axis=0)

perc =  array(perc).T

ld1 = 0.51374991
ld2 = 0.18869174

bfit = [ld1,ld2]
for p in perc:
  print p
  bfit += [p[2]]

print ' '
print 'best fit:'

print bfit

print ' '

triangle.corner(data)


vel1 = array(data_dict['eastvel'])
vel2 = array(data_dict['westvel'])

period = 2.2185733

RJ=69911.0

prad = 1.138*RJ

aprad = prad*(1.0 + atm)

apcircum = aprad*2*pi

rotvel = apcircum/(period*24*3600)

avel = (vel1 + vel2)/2

print percentile(avel,[0.15,2.275,15.865,50,84.135,97.725,99.85],axis=0)

dvel = (vel1 - vel2)

print percentile(dvel,[0.15,2.275,15.865,50,84.135,97.725,99.85],axis=0)

windvel = dvel - (rotvel*2)

print percentile(windvel,[0.15,2.275,15.865,50,84.135,97.725,99.85],axis=0)

savefig(openf.strip('.txt')+'.png')
savefig(openf.strip('.txt')+'.pdf')


#show()
hist(avel)
close

#show()
close()

hist(dvel)

#show()
close()

hist(windvel)

xlabel('Wind velocity')
ylabel('N')

title(openf.strip('.txt').replace('_',' '))

plt.gcf().subplots_adjust(bottom=0.15)

savefig(openf.strip('.txt')+'_windvel.png')
savefig(openf.strip('.txt')+'_windvel.pdf')

close()

hist(avel)

xlabel('Offset velocity')
ylabel('N')

title(openf.strip('.txt').replace('_',' '))

plt.gcf().subplots_adjust(bottom=0.15)

savefig(openf.strip('.txt')+'_offsetvel.png')
savefig(openf.strip('.txt')+'_offsetvel.pdf')

#show()


