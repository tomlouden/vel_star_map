# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

def add_diagram(ax,u1,u2):

  planet_radius = 0.1572
  b = 0.671

  atm = 0.3

  a_r = planet_radius*(1.0+planet_radius)

  gauplot([(0,0),(-b,-b),(0,-b),(b,-b)], [1.0,a_r,a_r,a_r], u1,u2,[-1,1], [-1,-0.33333])

  patches = []

  # Some limiting conditions on Wedge

  # sun patch
  #patches += [Wedge((0,0), 1.0, 0, 360,color=(1,1,0))]
  ingress = planet_patch(-b,-b,radius=planet_radius)
  patches += ingress

  full = planet_patch(0,-b,radius=planet_radius)
  patches += full

  egress = planet_patch(b,-b,radius=planet_radius)
  patches += egress

  for e in patches:
    ax.add_artist(e)

  ax.set_ylim(-1,-1+b)
  ax.set_xlim(-1,1)

  ax.get_xaxis().set_ticks([])
  ax.get_yaxis().set_ticks([])


def gauplot(centers, radiuses, u1,u2,xr=None, yr=None):
        nx, ny = 1000.,1000.
        xgrid, ygrid = np.mgrid[xr[0]:xr[1]:(xr[1]-xr[0])/nx,yr[0]:yr[1]:(yr[1]-yr[0])/ny]
        im = xgrid*0 + np.nan
        fis = np.concatenate((np.linspace(-np.pi,np.pi,100), [np.nan]) )
        thresh = 3
	i = 0
        for curcen,currad in zip(centers,radiuses):
                curim=(((xgrid-curcen[0])**2+(ygrid-curcen[1])**2)**.5)/currad*thresh
		if i == 0:
		  cmap = plt.cm.afmhot
		  cmap.set_bad('white')
		  color=(xgrid-curcen[0])[curim<thresh]
		  width = (abs(min(color.flatten()) - max(color.flatten())))*0.5
		  color = color/width
		  color = (color - min(color.flatten())) - 1.0
		  imcopy = im.copy()
		  imcopy[curim<thresh]=np.exp(-.5*curim**2)[curim<thresh]

		  r = curim[curim<thresh]
		  r = r / max(r.flatten())
		  theta = np.cos(np.arcsin(r))
		  intensity = (1-u1*(1-theta)-u2*(1-theta)**2)  # Fills the grid slice i at the x values generated by star with the star and the appropirate LD value.
		  imcopy[curim<thresh]=intensity
		  plt.imshow(imcopy.T, cmap=cmap, extent=xr+[-0.33333,-1])
		  i += 1
		else:
		  cmap = plt.cm.bwr
		  cmap.set_bad('white',alpha=0)
		  color=(xgrid-curcen[0])[curim<thresh]
		  width = (abs(min(color.flatten()) - max(color.flatten())))*0.5
		  color = color/width
		  color = (color - min(color.flatten())) - 1.0
		  im[curim<thresh]=np.sin(color*np.pi*0.5)
		  i += 1

        plt.imshow(im.T, cmap=cmap, extent=xr+[-0.33333,-1])

def planet_patch(x,y,radius=0.1,atm=1.1):
  patches = [
      Wedge((x,y), radius, 0, 360,color=(0,0,0)),             # Full circle
      #Wedge((x,y), radius*atm, -90, 90, width=radius*(atm-1),color=(1,1,1),alpha=1,linewidth=0), # Ring sector
      #Wedge((x,y), radius*atm, 90, -90, width=radius*(atm-1),color=(1,1,1),alpha=1,linewidth=0), # Ring sector
      #Wedge((x,y), radius*atm, -90, 90, width=radius*(atm-1),color=(1,0,0),alpha=0.5,linewidth=0), # Ring sector
      #Wedge((x,y), radius*atm, 90, -90, width=radius*(atm-1),color=(0,0,1),alpha=0.5,linewidth=0), # Ring sector
  ]
  return patches