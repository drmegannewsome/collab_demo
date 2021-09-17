#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 16:50:40 2021

@author: estefaniapadilla
"""

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import argparse
import numpy as np
import json
import pdb


#mpl.font_manager._rebuild()
mpl.rcParams['font.family'] = 'CMU Serif'
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.unicode_minus'] = False

markers = {'LCO':'o', 'ZTF':'*', 'Atlas':'s', 'Swift':'d'}


parser = argparse.ArgumentParser(description='plot LCO sn data prettily')
parser.add_argument('-f', '--filenames', type=str, help='name of data file', nargs='+', required=True)
parser.add_argument('-t', '--datatype', type=str, help='name of data type eg. LCO, ZTF, Swift, Atlas', nargs='+', required=True)
parser.add_argument('-o', '--output', type=str, help='name of output', required=False)
parser.add_argument('-x1', '--xmin', type=float, help='first date of plot', required=False)
parser.add_argument('-x2', '--xmax', type=float, help='last date of plot', required=False)
args = parser.parse_args()

# Put this in some visible spot where it can be easily added to
mjd_labels = {'mjd', 'bjd', 'jd'}
mag_labels = {'mag', 'magnitude', 'MAG'}
err_labels = {'err', 'mag_error', 'magnitude_error', 'error', 'dmag'}
filter_labels = {'filter', 'fltr', 'FILTER'}

#########


style = {'UVW2': {'color': '#FE0683', 'offset': -4.5},
         'UVM2': {'color': '#BF01BC', 'offset': -4.0},
         'UVW1': {'color': '#8B06FF', 'offset': -3.5},
         'U': {'color': '#3B0071', 'offset': -3.0},
	 'B': {'color': '#0057FF', 'offset': -3.0},
	 'g': {'color': '#00CCFF', 'offset': -2.0},
	 'V': {'color': '#78FF00', 'offset': -1.0},
	 'r': {'color': '#FF7C00', 'offset': 0.0},
	 'i': {'color': '#90002B', 'offset': +1.0},
         'c': {'color': 'magenta', 'offset': -1.5},
         'o': {'color': 'olive', 'offset': +2.0}}


alldata = {}
datatype = []
print(alldata)

for i in range(len(args.filenames)):
    if args.filenames[i][-4:] == 'json':
        tde_data = open(args.filenames[i])
        tde = json.load(tde_data)
        data = {'mjd': [], 'filter': [], 'mag': [], 'mag_error': []}
        for entry in tde['candidates']:
            if 'sigmapsf' in entry:
                data['mjd'].append(entry['mjd'])
                data['mag'].append(entry['magpsf'])
                data['mag_error'].append(entry['sigmapsf'])
                if entry['fid'] == 1:
                    data['filter'].append('r')
                elif entry['fid'] == 2:
                    data['filter'].append('g')
        table = Table(data)
        alldata[i] = table
        
    else:
    
        print(args.filenames[i])
        alldata[i] = ascii.read(args.filenames[i])
        datatype.append(args.datatype[i])
        
if args.output:
    output_name = args.output
else:
    output_name = 'temp_output'
if args.xmin and args.xmax:
    x1, x2 = args.xmin, args.xmax


fig, ax = plt.subplots()
if args.xmin != None :
    ax.set_xlim(x1, x2)
ax.invert_yaxis()



def get_labels(data): 
    # Find the intersection of the set of column names and the set of 
    # interchangeable labels
    #column_names = set(data.colnames)
    column_names= set(data.keys())
    mjd = column_names.intersection(mjd_labels).pop()
    mag = column_names.intersection(mag_labels).pop()
    err = column_names.intersection(err_labels).pop()
    filt = column_names.intersection(filter_labels).pop()
   
    return mjd, mag, err, filt, 




## Plot detections ##
def plot_detections(first_data, table, marker='o', alpha=0.75):
     for filters in style:
    
        data = first_data[(first_data[filt] == filters) &
                              (first_data[mag] < 100) & #good mags only
                              (first_data[err] < 100) ] #&
        if not data:
          continue
      
        if min(data[mjd]) > 2400000.5: # means it's a jd column
            date = data[mjd] -  2400000.5
            
        
        mjds, mags, errs = date, data[mag], data[err]
        for i in range(len(mjds)):
            table.add_row([mjds[i], mags[i], errs[i], filt, False])
        ax.errorbar(mjds, mags, yerr=errs,
                    alpha=alpha, marker=marker, markersize=8, linestyle='none', mec='black')


dtype = [('mjd', 'float'), ('mag','float'), ('err','float'), ('filter', 'str'), ('lim', 'bool')]
t = Table(data=np.zeros(1, dtype=dtype))

for i in range(len(args.filenames)):
    s =markers.get(datatype[i])
    mjd, mag, err, filt = get_labels(alldata[i])
    plot_detections(alldata[i], t , marker= s)
    

color_legend = []
tele_legend = []
for filt in style:
    if style[filt]['offset'] < 0:
        label = '{filt} {offset}'.format(filt=filt, offset=style[filt]['offset'])
    else:
        label = '{filt} +{offset}'.format(filt=filt, offset=style[filt]['offset'])
    color_legend.append(Line2D([0], [0], color=style[filt]['color'], label=label, lw=4))
for tele in markers:
    tele_legend.append(Line2D([0], [0], color='white', marker=markers[tele],
                              label=tele, mec='black', markersize=9))

ax.set_xlabel('MJD')
ax.set_ylabel('mag')
ax.set_title(output_name)
leg2 = ax.legend(handles=tele_legend, bbox_to_anchor=(1.01, 0.0), loc='lower left')
ax.add_artist(leg2)
ax.legend(handles=color_legend, bbox_to_anchor=(1.01, 1), loc='upper left')
plt.savefig('{name}.png'.format(name=output_name), dpi=300, bbox_inches='tight', overwrite=True)
plt.show()
plt.close()

ascii.write(t, 'alltelescopedata.txt', overwrite=True)
