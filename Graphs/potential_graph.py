import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, default='2DModel2')
args = parser.parse_args()

if args.model not in settings.MODELS:
    sys.exit(1)

cfg = settings.MODELS[args.model]
base_path = r"./MODEL"
output_path = r"./Graphs"
os.makedirs(output_path, exist_ok=True)

SMALL_SIZE = 14
MEDIUM_SIZE = 18
SIZE = 6

plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['axes.titlepad'] = 10 
plt.rcParams['axes.labelpad'] = 10 
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rc('font', size=SMALL_SIZE) 
plt.rc('axes', titlesize=MEDIUM_SIZE) 
plt.rc('axes', labelsize=MEDIUM_SIZE) 
plt.rc('xtick', labelsize=SMALL_SIZE) 
plt.rc('ytick', labelsize=SMALL_SIZE) 
plt.rc('legend', fontsize=12)

x_range = cfg.get('POT_x_range', [0.7, 10.5])
y_range = cfg.get('POT_y_range', [-4.0, 6.0])
nv_V0 = cfg.get('POT_nv_V0', 10)
nv_Vd = cfg.get('POT_nv_Vd', 10)
nv_VLCP = cfg.get('POT_nv_VLCP', 10)
n_ryd = cfg.get('POT_n_ryd', 5)

path_E_ryd = cfg.get('DR_Rydberg_E_path', '')
path_R_ryd = cfg.get('DR_Rydberg_R_path', '')

fig, ax = plt.subplots(figsize=(SIZE, SIZE))

def get_data(filename):
    path = os.path.join(base_path, filename)
    if os.path.exists(path):
        return np.loadtxt(path, comments='#')
    return None

def plot_levels(x_data, pot_data, energies, count, color):
    for i in range(min(count, len(energies))):
        e = energies[i]
        valid = pot_data <= e
        y_plot = np.full_like(x_data, np.nan)
        y_plot[valid] = e
        if i % 5 == 0:
            ax.plot(x_data, y_plot, color=color, linestyle='-', linewidth=0.9)
        else:
            ax.plot(x_data, y_plot, color=color, linestyle='--', dashes=(10, 5), linewidth=0.75)

data_V0 = get_data("V0.txt")
data_Vd = get_data("Vd.txt")
data_VLCP = get_data("VLCP.txt")

E_V0 = get_data("V0_eigE.txt")
E_Vd = get_data("Vd_eigE.txt")
E_VLCP = get_data("VLCP_eigE.txt")

if data_V0 is not None:
    x_V0, y_V0 = data_V0[:, 0], data_V0[:, 1]
    ax.plot(x_V0, y_V0, color='red', linestyle='-', linewidth=1.75, label=r'$V_0(R)$', zorder=5)
    if E_V0 is not None:
        plot_levels(x_V0, y_V0, E_V0[:, 1] if E_V0.ndim > 1 else E_V0, nv_V0, 'red')
        
    if os.path.exists(path_E_ryd):
        try:
            ryd_data = np.loadtxt(path_E_ryd)
            if os.path.exists(path_R_ryd):
                R_vals = np.loadtxt(path_R_ryd)
                if R_vals.ndim > 1:
                    R_vals = R_vals[:, 0]
                for i in range(min(n_ryd, ryd_data.shape[0])):
                    shift = np.interp(x_V0, R_vals, ryd_data[i, 1:])
                    ax.plot(x_V0, y_V0 + shift, color='gray', linestyle=':', linewidth=0.75, zorder=4)
            else:
                for i in range(min(n_ryd, ryd_data.shape[0])):
                    shift = ryd_data[i, 1] if ryd_data.ndim > 1 else ryd_data[i]
                    ax.plot(x_V0, y_V0 + shift, color='gray', linestyle=':', linewidth=0.75, zorder=4)
        except Exception:
            pass

if data_Vd is not None:
    x_Vd, y_Vd = data_Vd[:, 0], data_Vd[:, 1]
    ax.plot(x_Vd, y_Vd, color='green', linestyle='-', linewidth=1.75, label=r'$V_d(R)$', zorder=5)
    if E_Vd is not None:
        plot_levels(x_Vd, y_Vd, E_Vd[:, 1] if E_Vd.ndim > 1 else E_Vd, nv_Vd, 'green')

if data_VLCP is not None:
    x_VLCP, y_VLCP = data_VLCP[:, 0], data_VLCP[:, 1]
    gamma_VLCP = data_VLCP[:, 2]
    ax.plot(x_VLCP, y_VLCP, color='teal', linestyle='-.', linewidth=1.75, label=r'$V_\mathrm{LCP}(R)$', zorder=6)
    ax.fill_between(x_VLCP, y_VLCP + gamma_VLCP / 2, y_VLCP - gamma_VLCP / 2, color='green', alpha=0.4, zorder=1)
    if E_VLCP is not None:
        plot_levels(x_VLCP, y_VLCP, E_VLCP[:, 1] if E_VLCP.ndim > 1 else E_VLCP, nv_VLCP, 'teal')

ax.set_xlabel(r'Internuclear distance $R\,(a_0)$')
ax.set_ylabel(r'Potential energy $(\mathrm{eV})$')
ax.set_xlim(x_range)
ax.set_ylim(y_range)

ax.legend(frameon=False, loc='best')

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

fig.savefig(f'{output_path}/Potentials_{args.model}.pdf', format='pdf')
plt.close(fig)