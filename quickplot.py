# -*- coding: utf-8 -*-
from pylab import *
import triangle
import sys

nwalkers = 16

#openf = 'test_mcmc.txt'
#openf = 'mcmc_no_bootstrap.txt'

openf = sys.argv[-1]

#openf = 'correct_ld.txt'


data = []
i = 0
for line in open(openf):
  i += 1
  try:
    d = [float(x) for x in line.strip(' \n').split(' ')]
    if len(d) == 8:
      data += [d]
    else:
      print i, 'fail'
  except:
    print i, 'fail'

data = array(data)

print shape(data)

#data = loadtxt(openf)

print len(data), 'samples'
print len(data)/nwalkers, 'iterations'

burns = 4000

data = data[-1000:]

vels = data[:,2]

atm = data[:,4]

mask = [(vels > 2) & (vels < 20)]

data = data[mask]

perc = percentile(data,[2.275,15.865,50,84.135,97.725],axis=0)

perc =  array(perc).T

for p in perc:
  print p

triangle.corner(data)


vel1 = data[:,2]
vel2 = data[:,3]

period = 2.2185733

RJ=69911

prad = 1.138*RJ

aprad = prad*(1.0 + atm)

apcircum = aprad*2*pi

rotvel = apcircum/(period*24*3600)

print rotvel

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

#show()

