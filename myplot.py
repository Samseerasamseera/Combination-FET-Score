
import matplotlib.pyplot as plt
import numpy as np

def get_plot(pvl, cdf, k , type):
    plt.plot(pvl, cdf, marker='.', linestyle='none')

    x = np.arange(0, 1, 0.05)

    plt.xlabel('P-values')
    plt.set_xlim(x)
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Distribution Function (CDF) of P-values')
    plt.grid(True)
    plt.savefig(f'{k}_{type}_plot.svg',format='svg')


