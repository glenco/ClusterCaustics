#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 16:59:13 2018

@author: bmetcalf
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

def cum_dist(x,y = -1,inverse_order = False):
    '''
    Returns the cumulative distribution fit for 
    plotting.
    
    x : value for the x-axis, unordered
    y : value that is summed, ordered as x
        , if left out or -1 it will be the number of cases
    inverse_order : if False sum from the smallest 
        value of x, True from the largest
    '''
    if y is -1:
        y = np.ones(len(x))
    
    index = np.argsort(x)
 
    if(inverse_order) :
        index = index[::-1]
    
    y2 = y[index]
    x2 = x[index]
        
    y2 = np.repeat(np.cumsum(y2),2)
    x2 = np.roll(np.repeat(x2,2),-1)
    
    return x2[:-1],y2[:-1]
    
#df = pd.read_csv("DataFiles/snap_058_centered.txt.cy2049x2049.csv")
#df = pd.read_csv("DataFiles/snap_058_centered.txt.cy2049x2049S60.csv")

df_los = pd.read_csv("DataFiles/snap_058_centered.txt.cy2049x2049S60Zl0.506000LOSg.csv")

df_los = df_los[ df_los['caustic_area'] > 0]
df_los = df_los[ df_los['critical_area'] > 0]


df = pd.read_csv("DataFiles/snap_058_centered.txt.cy2049x2049S60Zl0.506000.csv")

df = df[ df['caustic_area'] > 0]
df = df[ df['critical_area'] > 0]


fact = (60*60*180/np.pi)**2
bins = 15

df1 = df[ df['type'] == 1] # radial
df2 = df[ df['type'] == 2] # tangential

x1 = fact * np.array( df1['caustic_area'])
x2 = fact * np.array( df2['caustic_area'])
x = fact * np.array( df['caustic_area'])

xx,y = cum_dist(x,-1,True)
#plt.fill_between(xx,y,label="radial")
plt.plot(xx,y,label="radial")

xx,y = cum_dist(x2,-1,True)
#plt.fill_between(xx,y,label="tangential")
plt.plot(xx,y,label="tangential")


df1_los = df_los[ df_los['type'] == 1] # radial
df2_los = df_los[ df_los['type'] == 2] # tangential

x1_los = fact * np.array( df1_los['caustic_area'])
x2_los = fact * np.array( df2_los['caustic_area'])
x_los = fact * np.array( df_los['caustic_area'])

xx,y = cum_dist(x_los,-1,True)
#plt.fill_between(xx,y,label="radial")
plt.plot(xx,y,label="with LOS")

xx,y = cum_dist(x2_los,-1,True)
plt.plot(xx,y,label="with LOS")


plt.xlim(1.0e-4,100)
plt.ylim(0,10)
plt.title("cumulative number of caustic curves")
plt.xscale('log')
plt.xlabel(r'area of caustic (arcsec$^2$)')
plt.ylabel(r'number of critical curves')
plt.legend()

#plt.savefig('caustic_area_distribution.png')
plt.show()

x1 = fact * np.array( df1['critical_area'] )
x2 = fact * np.array( df2['critical_area'] )
x = fact * np.array( df['critical_area'] )

xx,y = cum_dist(x,-1,True)
#plt.fill_between(xx,y,label="radial")
plt.plot(xx,y,label="radial")

xx,y = cum_dist(x2,-1,True)
#plt.fill_between(xx,y,label="tangential")
plt.plot(xx,y,label="tangential")


x1_los = fact * np.array( df1_los['critical_area'] )
x2_los = fact * np.array( df2_los['critical_area'] )
x_los = fact * np.array( df_los['critical_area'] )

xx,y = cum_dist(x_los,-1,True)
plt.plot(xx,y,label="with LOS")

xx,y = cum_dist(x2_los,-1,True)
plt.plot(xx,y,label="with LOS")

plt.xlim(1.0e-2,100)
plt.title("cumulative number of critical curves")
plt.xscale('log')
plt.xlabel(r'area of critical (arcsec$^2$)')
plt.ylabel(r'number of critical curves')
plt.legend()

#plt.savefig('critical_area_distribution.png')
plt.show()

#df1 = df1.sort_values(['critical_area'],ascending=0)
#df2 = df2.sort_values(['critical_area'],ascending=0)
#df = df.sort_values(['critical_area'],ascending=0)

#x1 = fact * df1['critical_area']
#x2 = fact * df2['critical_area']
#x = fact * df['critical_area']

y = fact * np.array(df['caustic_area'])

xx,yy = cum_dist(x,y,True)
#plt.fill_between(xx,yy,label="radial")
plt.plot(xx,yy,label="radial")

y2 = fact * np.array(df2['caustic_area'])

xx2,yy2 = cum_dist(x2,y2,True)
#plt.fill_between(xx2,yy2,label="tangential")
plt.plot(xx2,yy2,label="tangential")

y = fact * np.array(df_los['caustic_area'])
xx,yy = cum_dist(x_los,y,True)
plt.plot(xx,yy,label="with LOS")

y2 = fact * np.array(df2_los['caustic_area'])
xx2,yy2 = cum_dist(x2_los,y2,True)
plt.plot(xx2,yy2,label="with LOS")

plt.xscale('log')
plt.xlim(1.0e-2,4000)
#plt.ylim(143,160)
plt.ylim(260,285)
plt.ylim(240,310)

plt.ylabel(r'cumulative area within caustics (arcsec$^2$)')
plt.xlabel(r'area of critical (arcsec$^2$)')
plt.legend()

#plt.savefig('caust_area_tangential.png')
plt.show()

y = fact * np.array(df['caustic_area'])
xx,yy = cum_dist(x,y,True)
#plt.fill_between(xx,yy,label="radial")
plt.plot(xx,yy,label="radial")

y = fact * np.array(df_los['caustic_area'])
xx,yy = cum_dist(x_los,y,True)
plt.plot(xx,yy,label="with LOS")

plt.xscale('log')
plt.xlim(1.0e-2,4000)
plt.ylim(1550,1740)

plt.ylabel(r'cumulative area within caustics (arcsec$^2$)')
plt.xlabel(r'area of critical (arcsec$^2$)')
plt.legend()

#plt.savefig('caust_area_radial.png')
plt.show()

