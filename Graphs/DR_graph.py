import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator, LogLocator
import matplotlib.cm as cm
import gc
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings

parser = argparse.ArgumentParser(description="Script for generating DR cross section plots.")
parser.add_argument('--model', type=str, default='2DModel1', help='Select calculation mode')
args = parser.parse_args()

if args.model not in settings.MODELS:
    print(f"\n[ERROR] Mode '{args.model}' is not defined in settings.py!")
    sys.exit(1)

cfg = settings.MODELS[args.model]
print(f"--> Running {os.path.basename(__file__)} in mode: {args.model}")

base_path = r"./DATA"
output_path = r"./Graphs/DR"
os.makedirs(output_path, exist_ok=True)

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
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

path_NL = f"{base_path}/cs_DR.txt"
path_2D = cfg.get('DR_2D_path', '')
path_E_vib = r"./MODEL/V0_eigE.txt"
path_E_ryd = cfg.get('DR_Rydberg_E_path', '')
path_R_ryd = cfg.get('DR_Rydberg_R_path', '')

x_range = cfg.get('DR_x_range', [0, 1.5])
y_range = cfg.get('DR_y_range', [1e-6, 1e4])

x_range_c_list = cfg.get('DR_x_range_close', [[0.2, 0.4]])
if not isinstance(x_range_c_list[0], list): x_range_c_list = [x_range_c_list]

y_range_c_list = cfg.get('DR_y_range_close', [[1e-3, 1e3]])
if not isinstance(y_range_c_list[0], list): y_range_c_list = [y_range_c_list]

v_plotted_list = cfg.get('DR_v_plotted', [[1, 2, 3]])
if len(v_plotted_list) > 0 and not isinstance(v_plotted_list[0], list): 
    v_plotted_list = [v_plotted_list]

loc_main_list = cfg.get('DR_legend_loc', ['upper right'])
if isinstance(loc_main_list, str): loc_main_list = [loc_main_list]

loc_series_list = cfg.get('DR_legend_series_loc', ['upper left'])
if isinstance(loc_series_list, str): loc_series_list = [loc_series_list]

frame_main_list = cfg.get('DR_legend_frame', [True])
if isinstance(frame_main_list, bool): frame_main_list = [frame_main_list]

n_max = cfg.get('DR_n_max', 20)
annotate_n = cfg.get('DR_annotate_n', 4)
n_start = cfg.get('DR_n_start', 3)
leg_prefix = cfg.get('DR_legend_prefix', r'$\nu=')

cmap_name = cfg.get('DR_colormap', 'viridis')

def create_dr_plot(x_lim, y_lim, filename, plot_lines=False, loc_m='upper right', loc_s='upper left', frame_m=True, v_plot=[1,2,3]):
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
        data_NL = np.loadtxt(path_NL)
        x_NL, y_NL = data_NL[:, 0], data_NL[:, 1]
        mask = (x_NL >= x_lim[0]) & (x_NL <= x_lim[1]) & (y_NL > 0)
        line_nl, = ax.plot(x_NL[mask], y_NL[mask], color='red', linestyle='-', linewidth=1.0, label='Nonlocal', zorder=2)
        handles_main.append(line_nl)
        labels_main.append('Nonlocal')
    except OSError:
        pass

    try:
        data_2D = np.loadtxt(path_2D)
        x_2D, y_2D = data_2D[:, 0], data_2D[:, 1]
        mask = (x_2D >= x_lim[0]) & (x_2D <= x_lim[1]) & (y_2D > 0)
        line_2d, = ax.plot(x_2D[mask], y_2D[mask], color='blue', linestyle='-.', linewidth=0.7, label='Exact 2D', zorder=2)
        handles_main.append(line_2d)
        labels_main.append('Exact 2D')
    except OSError:
        pass

    lines_added = []
    labels_added = []
    x_margin = 0.05 * (x_lim[1] - x_lim[0]) 
    
    if plot_lines:
        try:
            vib_data = np.loadtxt(path_E_vib)
            E_vib = vib_data[:, 1]  
            E_v0 = E_vib[0]
            
            R_vals = np.loadtxt(path_R_ryd)
            if R_vals.ndim > 1: R_vals = R_vals[:, 0]
            r0_idx = (np.abs(R_vals - 2.0)).argmin()
            
            ryd_data = np.loadtxt(path_E_ryd)
            E_ryd = ryd_data[:, r0_idx + 1]

            valid_v_list = []
            peaks_data = {}

            for v in v_plot:
                if v >= len(E_vib): continue
                E_thresh = E_vib[v] - E_v0
                peaks = []
                
                if x_lim[0] <= E_thresh <= x_lim[1]:
                    peaks.append(('thresh', E_thresh))
                
                for idx in range(min(n_max, len(E_ryd))):
                    E_peak = E_thresh + E_ryd[idx]
                    if x_lim[0] <= E_peak <= x_lim[1]:
                        peaks.append(('ryd', E_peak, idx))
                
                if peaks:
                    valid_v_list.append(v)
                    peaks_data[v] = peaks

            if valid_v_list:
                cmap = plt.get_cmap(cmap_name)
                colors = [cmap(i) for i in np.linspace(0, 0.9, len(valid_v_list))]

                for i, v in enumerate(valid_v_list):
                    color = colors[i]
                    peaks = peaks_data[v]
                    
                    for peak_info in peaks:
                        if peak_info[0] == 'thresh':
                            ax.axvline(peak_info[1], color=color, linestyle='-', linewidth=1.2, alpha=0.8, zorder=1)
                        elif peak_info[0] == 'ryd':
                            E_peak = peak_info[1]
                            idx = peak_info[2]
                            ax.axvline(E_peak, color=color, linestyle=':', linewidth=0.8, alpha=0.8, zorder=1)
                            
                            if idx < annotate_n and E_peak > (x_lim[0] + x_margin):
                                n_label = n_start + idx
                                trans = ax.get_xaxis_transform()
                                ax.text(E_peak - 0.001*(x_lim[1]-x_lim[0]), 0.98, f"$n={n_label}$", 
                                        transform=trans, color=color, fontsize=8, 
                                        rotation=90, va='top', ha='right', zorder=3)

                    lines_added.append(plt.Line2D([], [], color=color, linestyle='-', linewidth=1.2))
                    labels_added.append(f'{leg_prefix}{v}$')

        except Exception as e:
            print(f"[WARNING] Could not plot structural lines: {e}")

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
    
    if plot_lines and lines_added:
        leg_series = ax.legend(lines_added, labels_added, loc=loc_s, 
                               frameon=frame_m, 
                               facecolor='white' if frame_m else 'none', 
                               edgecolor='black' if frame_m else 'none', 
                               framealpha=1.0 if frame_m else 0.0, 
                               fancybox=False, ncol=2, fontsize=10, handlelength=1, 
                               columnspacing=1.0, handletextpad=0.5)
        
        leg_series.get_frame().set_linewidth(0.5 if frame_m else 0.0)
        leg_series.set_zorder(10)
        
        if handles_main:
            ax.add_artist(leg_main)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
    ax.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

    fig.savefig(f'{output_path}/{filename}.pdf', format='pdf')
    fig.savefig(f'{output_path}/{filename}.pgf', format='pgf')
    plt.close(fig)
    gc.collect()

print("Generating Full Range DR plot...")
create_dr_plot(x_range, y_range, "DR_Full", plot_lines=False, loc_m='best', frame_m=False)

num_regions = len(x_range_c_list)
for i in range(num_regions):
    print(f"Generating Close-up DR plots (Region {i+1}/{num_regions})...")
    
    xc = x_range_c_list[i]
    yc = y_range_c_list[i] if i < len(y_range_c_list) else y_range_c_list[-1]
    vp = v_plotted_list[i] if i < len(v_plotted_list) else v_plotted_list[-1]
    lm = loc_main_list[i] if i < len(loc_main_list) else loc_main_list[-1]
    ls = loc_series_list[i] if i < len(loc_series_list) else loc_series_list[-1]
    fm = frame_main_list[i] if i < len(frame_main_list) else frame_main_list[-1]
    
    suffix = f"_{i+1}" if num_regions > 1 else ""

    create_dr_plot(xc, yc, f"DR_Close{suffix}", plot_lines=False, loc_m='best', frame_m=False)
    create_dr_plot(xc, yc, f"DR_Close_Lines{suffix}", plot_lines=True, loc_m=lm, loc_s=ls, frame_m=fm, v_plot=vp)

print(f"Done. Plots saved in {output_path}")