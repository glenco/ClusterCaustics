#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 13:00:19 2019

@author: bmetcalf
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pa

frac = False

filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj3def.csv'
name = 'def_prj3.png'
#filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj2def.csv'
#name = 'def_prj2.png'
#filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj1def.csv'
#name = 'def_prj1.png'

#filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj3fixdef.csv'
#name = 'fixdef_prj3.png'
#filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj2fixdef.csv'
#name = 'fixdef_prj2.png'
#filename = 'DataFiles/snap_058.sph1000x1000S30Zl0.506868Zs3.000000prj1fixdef.csv'
#name = 'fixdef_prj1.png'


#df = pa.read_csv(filename,sep=' ',names=['lens','image','x','y','dx','dy','mag'])
df = pa.read_csv(filename,sep='\s+')
#df = pa.read_csv(filename,sep=' ',nrows=65)

for col in df.columns :
    print(col)


radToArcs = 180.*60*60/np.pi

df['delta_x'] = radToArcs*df['delta_x']
df['delta_y'] = radToArcs*df['delta_y']
df['x_image'] = radToArcs*df['x_image']
df['y_image'] = radToArcs*df['y_image']

#plt.scatter(df['delta_x'],df['delta_y'],s=0.7,alpha=0.)a

neg_images = df[ df['mag'] < 0]
plt.scatter(neg_images['delta_x'],neg_images['delta_y'],s=0.7,alpha=1.0)

pos_images = df[ df['mag'] > 0]
plt.scatter(pos_images['delta_x'],pos_images['delta_y'],s=0.7,alpha=0.7)

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
    same_number = data['same_number'].iloc[0]
    #change_image_number = ( data.at[0,'same_number'] == 0 )
    
    if( n > 1 and same_number ) :
        xs = np.array(data['x_image'])
        ys = np.array(data['y_image'])

        dxs = np.array(data['delta_x'])
        dys = np.array(data['delta_y'])

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
    
                xx[ii] = dxn/s
                yy[ii] = dyn/s
                #xx[ii] = abs(dxn)
                #yy[ii] = dyn
            
                ii += 1

xx = xx[0:ii]
yy = yy[0:ii]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(xx,yy,s=1.8,alpha=0.7)

if frac :
    plt.xlim(-0.5,0.5)
    plt.ylim(0,0.5)
    
    plt.xlabel(r'$\delta_\parallel / \Delta \theta$')
    plt.ylabel(r'$\delta_\perp /  \Delta \theta$')

else :
    plt.xlim(-0.5,0.5)
    plt.ylim(0,0.5)
 
    plt.xlabel('arcsec')
    plt.ylabel('arcsec')
    
ax.set_aspect(aspect=1.0)

s = xx * xx + yy * yy

s = np.sort(s)

r = np.sqrt( s[ int(0.95*ii) ] )
x = r * np.cos( np.arange(0,np.pi,np.pi/1000) )
y = r * np.sin( np.arange(0,np.pi,np.pi/1000) )

plt.plot(x,y,label='95%')

r = np.sqrt( s[ int(0.68*ii) ] )
x = r * np.cos( np.arange(0,np.pi,np.pi/1000) )
y = r * np.sin( np.arange(0,np.pi,np.pi/1000) )

plt.plot(x,y,label='68%')

plt.legend()
plt.savefig(name)
plt.show()