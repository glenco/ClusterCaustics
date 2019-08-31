#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 13:00:19 2019

@author: bmetcalf
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pa

#filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj3def.csv'
filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj1def.csv'

df = pa.read_csv(filename,sep=' ',names=['lens','image','x','y','dx','dy','mag'])

radToArcs = 180.*60*60/np.pi

df['dx'] = radToArcs*df['dx']
df['dy'] = radToArcs*df['dy']
df['x'] = radToArcs*df['x']
df['y'] = radToArcs*df['y']

#plt.scatter(df['dx'],df['dy'],s=0.7,alpha=0.)

neg_images = df[ df['mag'] < 0]
plt.scatter(neg_images['dx'],neg_images['dy'],s=0.7,alpha=1.0)

pos_images = df[ df['mag'] > 0]
plt.scatter(pos_images['dx'],pos_images['dy'],s=0.7,alpha=0.7)

plt.xlim(-10,10)
plt.ylim(-10,10)

plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.show()


Nlens = df.at[df.shape[0] - 1,'lens'] + 1

xx = np.zeros(Nlens*6 )
yy = np.zeros(Nlens*6 )
ii = 0
for l in range(Nlens) :
    data = df[ df['lens'] == l]
    n = data.shape[0]

    xs = np.array(data['x'])
    ys = np.array(data['y'])

    dxs = np.array(data['dx'])
    dys = np.array(data['dy'])

    if n > 1 :
        for i in range(n) :
            for j in range(i+1,n) :
                x = xs[i] - ys[j]
                y = ys[i] - ys[j]
                
                s = np.sqrt(x*x + y*y)
                x = x/s
                y = y/s
    
                dx = dxs[i] - dxs[j]
                dy = dys[i] - dys[j]
                
                ds = np.sqrt(dx*dx + dy*dy)
   
                # deflection in frame of seporation
                dxn = x*dx + y*dy
                dyn = np.sqrt( ds*ds  - dxn*dxn )    
    
                xx[ii] = abs(dxn)/s
                yy[ii] = dyn/s
                #xx[ii] = abs(dxn)
                #yy[ii] = dyn
            
                ii += 1

xx = xx[0:ii]
yy = yy[0:ii]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(xx,yy,s=1.8,alpha=0.7)

plt.xlim(0,0.5)
plt.ylim(0,0.5)
plt.xlabel('arcsec')
plt.ylabel('arcsec')
ax.set_aspect(aspect=1.0)

s = xx * xx + yy * yy
print s

s = np.sort(s)

r = np.sqrt( s[ int(0.95*ii) ] )
x = r * np.cos( np.arange(0,np.pi/2,np.pi/1000) )
y = r * np.sin( np.arange(0,np.pi/2,np.pi/1000) )

plt.plot(x,y,label='95%')

r = np.sqrt( s[ int(0.68*ii) ] )
x = r * np.cos( np.arange(0,np.pi/2,np.pi/1000) )
y = r * np.sin( np.arange(0,np.pi/2,np.pi/1000) )

plt.plot(x,y,label='68%')

plt.legend()
plt.show()