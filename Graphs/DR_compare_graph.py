import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator, LogLocator
import gc
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings

parser = argparse.ArgumentParser(description="Script for generating DR cross section comparison plots.")
parser.add_argument('--model', type=str, default='2DModel2')
args = parser.parse_args()

if args.model not in settings.MODELS:
    sys.exit(1)

cfg = settings.MODELS[args.model]

base_path = r"./DATA"
output_path = r"./Graphs/DRcomparison"
os.makedirs(output_path, exist_ok=True)

SMALL_SIZE = 14
MEDIUM_SIZE = 18
SIZE = 6
aspect_ratio = 'golden_ratio'

plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['axes.titlepad'] = 10 
plt.rcParams['axes.labelpad'] = 10 
plt.rcParams["figure.autolayout"] = True   
plt.rc('font', size=SMALL_SIZE) 
plt.rc('axes', titlesize=MEDIUM_SIZE) 
plt.rc('axes', labelsize=MEDIUM_SIZE) 
plt.rc('xtick', labelsize=SMALL_SIZE) 
plt.rc('ytick', labelsize=SMALL_SIZE) 

X_axis = r'Incoming electron energy $\epsilon_i\,(\mathrm{eV})$'
Y_axis = r'DR cross section $\sigma_\mathrm{DR}\,(\mathrm{a.u.})$'

path_data = f"{base_path}/cs_DR.txt"

x_range = cfg.get('DR_x_range', [0, 1.5])
y_range = cfg.get('DR_y_range', [1e-6, 1e4])

x_range_c_list = cfg.get('DR_x_range_close', [[0.2, 0.4]])
if not isinstance(x_range_c_list[0], list): x_range_c_list = [x_range_c_list]

y_range_c_list = cfg.get('DR_y_range_close', [[1e-3, 1e3]])
if not isinstance(y_range_c_list[0], list): y_range_c_list = [y_range_c_list]

loc_main_list = cfg.get('DR_legend_loc', ['upper right'])
if isinstance(loc_main_list, str): loc_main_list = [loc_main_list]

frame_main_list = cfg.get('DR_legend_frame', [True])
if isinstance(frame_main_list, bool): frame_main_list = [frame_main_list]

def create_drc_plot(x_lim, y_lim, filename, loc_m='upper right', frame_m=True):
    if aspect_ratio == 'golden_ratio':
        golden = (1 + 5 ** 0.5) / 2
        fig, ax = plt.subplots(figsize=(SIZE * golden, SIZE))
    else:
        fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    ax.set_yscale('log')
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    handles_main = []
    labels_main = []
    
    try:
        data = np.loadtxt(path_data)
        x = data[:, 0]
        y1 = data[:, 1]
        y2 = data[:, 2]
        y3 = data[:, 3]
        
        mask1 = (x >= x_lim[0]) & (x <= x_lim[1]) & (y1 > 0)
        line1, = ax.plot(x[mask1], y1[mask1], color='red', linestyle='-', linewidth=0.8, label='T-matrix', zorder=2)
        handles_main.append(line1)
        labels_main.append('T-matrix')
        
        mask2 = (x >= x_lim[0]) & (x <= x_lim[1]) & (y2 > 0)
        line2, = ax.plot(x[mask2], y2[mask2], color='blue', linestyle='-.', linewidth=0.6, label='Outgoing flux', zorder=3)
        handles_main.append(line2)
        labels_main.append('Outgoing flux')
        
        mask3 = (x >= x_lim[0]) & (x <= x_lim[1]) & (y3 > 0)
        line3, = ax.plot(x[mask3], y3[mask3], color='green', linestyle='--', linewidth=0.4, label='Optical theorem', zorder=4)
        handles_main.append(line3)
        labels_main.append('Optical theorem')
        
    except OSError:
        pass

    ax.set_xlabel(X_axis)
    ax.set_ylabel(Y_axis)
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    
    if handles_main:
        leg_main = ax.legend(handles_main, labels_main, loc=loc_m, 
                             frameon=frame_m, 
                             facecolor='white' if frame_m else 'none', 
                             edgecolor='black' if frame_m else 'none', 
                             framealpha=1.0 if frame_m else 0.0, 
                             fancybox=False, ncol=1, fontsize=12, handlelength=2)
        
        leg_main.get_frame().set_linewidth(0.5 if frame_m else 0.0)
        leg_main.set_zorder(10)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
    ax.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

    fig.savefig(f'{output_path}/{filename}.pdf', format='pdf')
    fig.savefig(f'{output_path}/{filename}.pgf', format='pgf')
    plt.close(fig)
    gc.collect()

create_drc_plot(x_range, y_range, "DRC_Full", loc_m='best', frame_m=False)

num_regions = len(x_range_c_list)
for i in range(num_regions):
    xc = x_range_c_list[i]
    yc = y_range_c_list[i] if i < len(y_range_c_list) else y_range_c_list[-1]
    lm = loc_main_list[i] if i < len(loc_main_list) else loc_main_list[-1]
    fm = frame_main_list[i] if i < len(frame_main_list) else frame_main_list[-1]
    
    suffix = f"_{i+1}" if num_regions > 1 else ""

    create_drc_plot(xc, yc, f"DRC_Close{suffix}", loc_m=lm, frame_m=fm)