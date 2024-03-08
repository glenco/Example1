#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 11:53:12 2024

This is for plotting caustics and critical curves for GLAMER example 1

@author: bmetcalf
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pa

df = pa.read_csv('caustics_and_crits.csv',comment='#')

ids = set(df['id'])
print(ids)

colors = ['red','blue','green']

for i in ids :
    
    df2 = df[df['id'] == i]
    plt.plot(df2['crit_x'],df2['crit_y'],color=colors[i])
    plt.plot(df2['caust_x'],df2['caust_y'],color=colors[i])
    
plt.show()