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
    
df = pd.read_csv("DataFiles/snap_058_centered.txt.cy2049x2049.csv")

df = df[ df['caustic_area'] > 0]
df = df[ df['critical_area'] > 0]

df1 = df[ df['type'] == 1] # radial
df2 = df[ df['type'] == 2] # tangential

fact = (60*60*180/np.pi)**2
bins = 15

x1 = fact * np.array( df1['caustic_area'])
x2 = fact * np.array( df2['caustic_area'])
x = fact * np.array( df['caustic_area'])

xx,y = cum_dist(x,-1,True)
plt.fill_between(xx,y,label="radial")

xx,y = cum_dist(x2,-1,True)
plt.fill_between(xx,y,label="tangential")

plt.xlim(1.0e-4,100)
plt.title("cumulative number of caustic curves")
plt.xscale('log')
plt.xlabel(r'area of caustic (arcsec$^2$)')
plt.ylabel(r'number of critical curves')
plt.legend()

plt.savefig('caustic_area_distribution.png')
plt.show()

x1 = fact * np.array( df1['critical_area'] )
x2 = fact * np.array( df2['critical_area'] )
x = fact * np.array( df['critical_area'] )

xx,y = cum_dist(x,-1,True)
plt.fill_between(xx,y,label="radial")

xx,y = cum_dist(x2,-1,True)
plt.fill_between(xx,y,label="tangential")

plt.xlim(1.0e-2,100)
plt.title("cumulative number of critical curves")
plt.xscale('log')
plt.xlabel(r'area of critical (arcsec$^2$)')
plt.ylabel(r'number of critical curves')
plt.legend()

plt.savefig('critical_area_distribution.png')
plt.show()

#df1 = df1.sort_values(['critical_area'],ascending=0)
#df2 = df2.sort_values(['critical_area'],ascending=0)
#df = df.sort_values(['critical_area'],ascending=0)

#x1 = fact * df1['critical_area']
#x2 = fact * df2['critical_area']
#x = fact * df['critical_area']

y = fact * np.array(df['caustic_area'])

xx,yy = cum_dist(x,y,True)
plt.fill_between(xx,yy,label="radial")

y2 = fact * np.array(df2['caustic_area'])

xx2,yy2 = cum_dist(x2,y2,True)
plt.fill_between(xx2,yy2,label="tangential")

plt.xscale('log')
plt.xlim(1.0e-2,4000)
plt.ylim(143,160)

plt.ylabel(r'cumulative area within caustics (arcsec$^2$)')
plt.xlabel(r'area of critical (arcsec$^2$)')
plt.legend()

plt.savefig('caust_area_tangential.png')
plt.show()

plt.fill_between(xx,yy,label="radial")

plt.xscale('log')
plt.xlim(1.0e-2,4000)
plt.ylim(945,1030)

plt.ylabel(r'cumulative area within caustics (arcsec$^2$)')
plt.xlabel(r'area of critical (arcsec$^2$)')
plt.legend()

plt.savefig('caust_area_radial.png')
plt.show()

