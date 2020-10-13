import matplotlib.pyplot as plt
import numpy as np


def plot_x_y(ax, x_values, y_values, x_name, y_name, line_type,
             **kwargs):
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)

    ax.ticklabel_format(axis="x", style="plain")
    ax.plot(x_values, y_values, line_type, **kwargs)



