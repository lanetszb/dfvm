import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from utilities import plot_x_y
from matplotlib import rc
import matplotlib.lines as mlines

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

# Figure validating the diffusion numerical solution against comsole
table_coms = pd.read_csv('../inOut/conc_coord_comsole.csv')
df_coms = pd.DataFrame(table_coms)
#
table_dfvm = pd.read_csv('../inOut/conc_coord_dfvm.csv')
df_dfvm = pd.DataFrame(table_dfvm)
#
fig_width = 4.5
y_scale = 0.9
#

fig1, ax1 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)
#
coms_arr = [np.array(df_coms.coord_coms),
            np.array(df_coms.conc_coms1),
            np.array(df_coms.conc_coms2),
            np.array(df_coms.conc_coms3),
            np.array(df_coms.conc_coms4)]

dfvm_arr = [np.array(df_dfvm.coord_dfvm),
            np.array(df_dfvm.conc_dfvm1),
            np.array(df_dfvm.conc_dfvm2),
            np.array(df_dfvm.conc_dfvm3),
            np.array(df_dfvm.conc_dfvm4)]
#
colors = ['b', 'r', 'g', 'y']

for i in range(len(coms_arr) - 1):
    plot_x_y(ax1, coms_arr[0], coms_arr[i + 1], x_name='Length, m',
             y_name='Concentration, $kg/m^3$', line_type='-',
             color=colors[i],
             lw=1.5)

for i in range((len(dfvm_arr) - 1)):
    plot_x_y(ax1, dfvm_arr[0], dfvm_arr[i + 1], x_name='Length, $m$',
             y_name='Concentration, $kg/m^3$',
             line_type='--',
             dashes=(5, 10),
             lw=1.5,
             color='k',
             marker='',
             markersize=1.5)

markers = ['', '']
lws = [1.5, 1.5]
line_names = ['solver', 'comsol']
line_types = ['--', '-']
for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], marker=markers[i], lw=lws[i],
             markersize='2.5')

legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc=1)

plt.gca().add_artist(legend_1)

ax2 = mlines.Line2D([], [], linestyle='-', c='b', markersize='2.0',
                    label=str(0.01))
ax3 = mlines.Line2D([], [], linestyle='-', c='r', markersize='2.0',
                    label=str(0.06))
ax4 = mlines.Line2D([], [], linestyle='-', c='g', markersize='2.0',
                    label=str(0.25))
ax5 = mlines.Line2D([], [], linestyle='-', c='g', markersize='2.0',
                    label=str(0.75))

plt.legend(handles=[ax2, ax3, ax4, ax5], scatterpoints=1, frameon=True,
           labelspacing=0.5,
           title='time, sec',
           loc=9)

plt.show()
#
plt.savefig('../inOut/validation_comsole.eps', format="eps",
            bbox_inches='tight')

# Figure validating the diffusion numerical and analytical solutions
table_num = pd.read_csv('../inOut/conc_time_num.csv')
df_num = pd.DataFrame(table_num)

table_analyt = pd.read_csv('../inOut/conc_time_analyt.csv')
df_analyt = pd.DataFrame(table_analyt)

fig_width = 4.5
y_scale = 0.9
#

fig2, ax2 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)

num_arr = [np.array(df_num.coord_num),
           np.array(df_num.conc_num1),
           np.array(df_num.conc_num2),
           np.array(df_num.conc_num3),
           np.array(df_num.conc_num4)]

analyt_arr = [np.array(df_analyt.coord_analyt),
              np.array(df_analyt.conc_analyt1),
              np.array(df_analyt.conc_analyt2),
              np.array(df_analyt.conc_analyt3),
              np.array(df_analyt.conc_analyt4)]

colors = ['b', 'r', 'g', 'y']

for i in range(len(analyt_arr) - 1):
    plot_x_y(ax2, analyt_arr[0], analyt_arr[i + 1], x_name='Length, m',
             y_name='Concentration, $kg/m^3$',
             line_type='-',
             color=colors[i],
             lw=1.5)

for i in range(len(num_arr) - 1):
    plot_x_y(ax2, num_arr[0], num_arr[i + 1], x_name='Length, $m$',
             y_name='Concentration, $kg/m^3$',
             line_type='--',
             dashes=(5, 10),
             lw=1.5,
             color='k',
             marker='',
             markersize=1.5)

markers = ['', '']
lws = [1.5, 1.5]
line_names = ['numerical', 'analytical']
line_types = ['--', '-']
for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], marker=markers[i], lw=lws[i],
             markersize='2.5')

legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc=1)

plt.gca().add_artist(legend_1)

ax3 = mlines.Line2D([], [], linestyle='-', c='b', markersize='2.0',
                    label=str(0.01))
ax4 = mlines.Line2D([], [], linestyle='-', c='r', markersize='2.0',
                    label=str(0.03))
ax5 = mlines.Line2D([], [], linestyle='-', c='g', markersize='2.0',
                    label=str(0.08))
ax6 = mlines.Line2D([], [], linestyle='-', c='g', markersize='2.0',
                    label=str(0.18))

plt.legend(handles=[ax3, ax4, ax5, ax6], scatterpoints=1, frameon=True,
           labelspacing=0.5,
           title='time, sec',
           loc=9)

plt.savefig('../inOut//validation_analytics.eps', format="eps",
            bbox_inches='tight')
