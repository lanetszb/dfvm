import matplotlib.pyplot as plt
import numpy as np


def plot_x_y(x_values, y_values, x_name, y_name, graph_name, line_type,
             **kwargs):
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    plt.suptitle(graph_name)

    plt.xlabel(x_name)
    plt.ylabel(y_name)

    plt.ticklabel_format(axis="x", style="plain")
    # ax = plt.gca()
    plt.plot(x_values, y_values, line_type, **kwargs)
    plt.legend(y_name)
    plt.show()
    # ax.yaxis.set_major_locator(plt.LinearLocator(numticks=6))


